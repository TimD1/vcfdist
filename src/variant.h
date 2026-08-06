/**
 * @file variant.h
 * @brief Per-contig and per-callset variant containers with VCF parsing and output utilities.
 */
#ifndef _VARIANT_H_
#define _VARIANT_H_

#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#include "htslib/vcf.h"

#include "fasta.h"
#include "defs.h"

/**
 * @struct hap_fields
 * @brief One haplotype's evaluation results for a single variant, set during prec_recall_aln().
 */
struct hap_fields {
    errtype_t errtype = ERRTYPE_UN; ///< error type: TP, FP, FN
    int sync_group = 0;           ///< group of variants that participate in credit
    float callq = 0;              ///< min call quality in sync group
    int ref_ed = 0;               ///< reference edit distance in sync group
    int query_ed = 0;             ///< query edit distance in sync group
    float credit = 0;             ///< percentage reduction in edit dist (ref->query)
};

/**
 * @struct var_fields
 * @brief Every field describing a single variant, named rather than positional.
 *
 * Members are matched by designator at each call site, so inserting a field cannot rebind the
 * others. The members with no default initializer are required: omitting one is a build failure
 * under -Werror=missing-field-initializers, while omitting a defaulted member is not. Designators
 * must appear in declaration order, which follows the storage order in ctgVariants.
 *
 * Call sites name the type -- add_var(var_fields{.pos = ...}), not add_var({.pos = ...}) -- because
 * GCC 13.3 rejects a bare designated-initializer list as a function argument once any member is
 * initialized from a non-constant expression. Naming the type still enforces the required members.
 */
struct var_fields {
    int pos;                      ///< variant start position (0-based)
    int rlen;                     ///< reference length
    edittype_t type;              ///< variant type: NONE, SUB, INS, DEL, CPX
    bedloc_t loc;                 ///< BED location: INSIDE, OUTSIDE, BORDER
    std::string ref;              ///< variant reference allele
    std::string alt;              ///< variant alternate allele
    gt_t orig_gt;                 ///< simple genotype (0|1, 1|0, or 1|1)
    float gt_qual;                ///< genotype quality (capped above at --max-qual when stored)
    float var_qual;               ///< variant quality (capped above at --max-qual when stored)
    int phase_set;                ///< integer representing variant phase set (0 = missing)
    int rec_idx = -1;             ///< source VCF record ordinal (0-based, -1 = unknown)
    int alt_idx = -1;             ///< original ALT ordinal (1-based, -1 = unknown)
    ploidy_t ploidy = PLOIDY_DIPLOID; ///< variant ploidy, from std::abs(ngt)
    int supercluster = -1;        ///< supercluster index (-1 = not yet assigned)
    gt_t matched_gt = GT_REF_REF; ///< the other callset's genotype, recovered by alignment
    EnumArray<hap_t, hap_fields, HAP_SLOTS> hap = {}; ///< per-haplotype results, indexed by HAP1 and HAP2
};

/**
 * @class srcRecords
 * @brief Site-level columns retained from one input VCF, indexed by source record ordinal.
 *
 * The summary VCF is synthesized from the internal variant arrays rather than copied from the
 * input, so the columns it cannot derive are held here for the run's lifetime and read back by
 * the writer. Indexing by record ordinal rather than by variant stores one copy per source
 * record, which the two halves of a split complex variant and the two entries of a het-alt
 * record share; the entries of a het-alt record then subset that copy's ALT-indexed fields to
 * their own allele. Columns are retained before any filtering decision, since a record excluded
 * from evaluation is still written out, so every vector is indexable by any ordinal below its size.
 *
 * Columns are held as the VCF text htslib rendered them from the input, and are re-typed against
 * the output header's own declarations when written. Holding one string per column per record
 * costs far less than the tagged per-field values a typed store would need, and both directions
 * go through htslib rather than hand-rolled formatting.
 */
class srcRecords {
public:

    /** @brief Stores one record's preserved columns at its 0-based ordinal within the input VCF. */
    void add(int rec_idx, const std::string & id, const std::string & qual,
            const std::string & filter, const std::string & info,
            const std::string & fmt_keys, const std::string & fmt_vals);

    /** @brief Releases the spare capacity that appending record by record left behind. */
    void shrink();

    std::vector<std::string> ids;       ///< ID column, verbatim
    std::vector<std::string> quals;     ///< QUAL column, verbatim
    std::vector<std::string> filters;   ///< FILTER column, verbatim
    std::vector<std::string> infos;     ///< INFO column, preserved fields only
    std::vector<std::string> fmt_keys;  ///< preserved FORMAT keys, each prefixed with ':'
    std::vector<std::string> fmt_vals;  ///< sample values parallel to fmt_keys, each ':'-prefixed

    // set once from the input VCF header (size equal to each other)
    std::vector<std::string> hdr_keys;  ///< "<line type>/<ID>" of each retained header line
    std::vector<std::string> hdr_lines; ///< FILTER, INFO, and FORMAT lines declaring those fields

    // set once from the input VCF header; a key absent from either map is not ALT-indexed
    std::unordered_map<std::string, int> ///< BCF_VL_A/R/G length class of each INFO key
        info_lens;
    std::unordered_map<std::string, int> ///< BCF_VL_A/R/G length class of each FORMAT key
        fmt_lens;
};

/**
 * @struct sample_fields
 * @brief One sample's FORMAT values for a summary VCF record, in FORMAT declaration order.
 *
 * Values are stored htslib-encoded: a per-allele entry that was never evaluated holds
 * bcf_int32_missing or a bcf_float_set_missing() float, and the two Number=. String fields hold the
 * comma-joined list htslib stores as a single string. A sample that made no call at this locus
 * leaves the per-allele vectors empty, which writes '.' for every one of its fields.
 *
 * The fields carried over from the input follow the fixed ones, and only the sample owning the
 * record supplies them; the other sample's stay empty and are written as missing.
 */
struct sample_fields {
    std::vector<int32_t> gt;   ///< GT alleles, bcf_gt_phased()-encoded, in GT allele order
    std::string bd = ".";      ///< BD, per-allele benchmark decision (TP/FP/FN)
    std::vector<float> bc;     ///< BC, per-allele benchmark credit
    std::vector<int32_t> rd;   ///< RD, per-allele reference edit distance
    std::vector<int32_t> qd;   ///< QD, per-allele query edit distance
    std::string bk = ".";      ///< BK, per-allele benchmark category ('gm', 'lm', or '.')
    float qq = 0;              ///< QQ, variant quality
    int32_t sc = 0;            ///< SC, supercluster index in contig
    std::vector<int32_t> sg;   ///< SG, per-allele sync group
    int32_t ps = 0;            ///< PS, input phase set
    int32_t pb = 0;            ///< PB, output phase block
    int32_t bs = 0;            ///< BS, block phase
    int32_t vp = 0;            ///< VP, variant phase
    int32_t fe = 0;            ///< FE, flip error
    std::string ge = ".";      ///< GE, allele count (genotype) error
    std::string src_keys;      ///< preserved source FORMAT keys, each prefixed with ':'
    std::string src_vals;      ///< values parallel to src_keys, each prefixed with ':'
};

/**
 * @class ctgSideline
 * @brief Variants of one contig and callset retained in the output but excluded from evaluation.
 *
 * A retained variant never enters ctgVariants, so every loop over `poss`/`n` in clustering,
 * alignment, phasing, and reporting is unreachable for it and "excluded from all analysis" holds by
 * construction rather than by a check at each of those sites. The summary VCF writer is the only
 * consumer: it merges these entries into its position-ordered walk and writes each as a record no
 * callset was evaluated on.
 *
 * Entries are keyed by (source record ordinal, haplotype), since a per-allele reason can leave one
 * haplotype evaluated while the other is not. A record-scope reason is decided before the genotype
 * is read, so it applies to the whole record and stores SIDELINE_ALL_HAPS.
 *
 * The site columns come from the source record verbatim, because a retained record was never
 * normalized or split: nothing derived them the way it does for an evaluated variant. Its whole ALT
 * list is emitted as written, so an ALT-indexed value still indexes the alleles beside it and is
 * carried over unsubset, at the cardinality the propagated declaration already states. Nothing here
 * subsets or reformats a retained column, so the store is read where it is written rather than
 * through the per-column accessors ctgVariants needs for its own subsetting.
 */
class ctgSideline {
public:

    /** @brief Constructs a contig-specific container of retained variants. */
    ctgSideline(const std::string & ctg);

    /** @brief Appends one retained variant, in source record order. */
    void add(int rec_idx, int8_t hap, int pos, const std::string & ref, const std::string & alt,
            const std::string & gt, sideline_t reason);

    /** @brief Sets the fixed VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) of one record. */
    void set_var_record(const bcf_hdr_t* hdr, bcf1_t* rec, const std::string & ctg, int si) const;

    /** @brief Returns the FORMAT values of the sample that called one retained variant. */
    sample_fields var_sample_fields(int si) const;

    std::string ctg;                 ///< Contig name (chromosome identifier)
    std::vector<int> rec_idxs;       ///< source VCF record ordinal (0-based)
    std::vector<int8_t> haps;        ///< haplotype the reason applies to (SIDELINE_ALL_HAPS = record)
    std::vector<int> poss;           ///< source record start position (0-based), for output order
    std::vector<std::string> refs;   ///< REF column, verbatim
    std::vector<std::string> alts;   ///< ALT column, verbatim (comma-separated if multi-allelic)
    std::vector<std::string> gts;    ///< sample's GT value, verbatim ("." if the record has none)
    std::vector<sideline_t> reasons; ///< retention reason, selecting this record's FILTER tag
    int n = 0;                       ///< Total number of retained variants

    // shared with every other container of this callset, indexed by rec_idxs
    std::shared_ptr<srcRecords> src_recs; ///< retained source records (nullptr = none retained)
};

/**
 * @class ctgVariants
 * @brief Store all variant information for a single contig and callset.
 */
class ctgVariants {
public:

    /** @brief Constructs a contig-specific variant container with empty data structures. */
    ctgVariants(const std::string & ctg);

    /** @brief Appends a variant with all fields explicitly specified. */
    void add_var(const var_fields & var);

    /** @brief Returns every field of one variant, for copying it into another container. */
    var_fields get_var(int idx) const;

    /** @brief Sets the fixed VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER) of one record. */
    void set_var_record(const bcf_hdr_t* hdr, bcf1_t* rec, std::shared_ptr<fastaData> ref,
            const std::string & ctg, int idx) const;

    /** @brief Returns one sample's FORMAT values for a variant it called. */
    sample_fields var_sample_fields(int vi, int sc_idx, int phase_block,
            bool phase_switch, bool phase_flip, bool query = false,
            bool owns_record = false) const;

    /** @brief Returns the source record's ID column, or "." if none was retained. */
    const std::string & src_id(int vi) const;

    /** @brief Returns the source record's QUAL column, or "." if none was retained. */
    const std::string & src_qual(int vi) const;

    /** @brief Returns the source record's FILTER column, or "PASS" if none was retained. */
    const std::string & src_filter(int vi) const;

    /** @brief Returns the source record's INFO column subset to this variant's allele, or ".". */
    std::string src_info(int vi) const;

    /** @brief Returns the source record's preserved FORMAT keys, each prefixed with ':'. */
    const std::string & src_fmt_keys(int vi) const;

    /** @brief Returns the source sample's FORMAT values subset to this variant's allele. */
    std::string src_fmt_vals(int vi) const;

    /** @brief Returns true if a variant is present on the specified haplotype. */
    bool var_on_hap(int var_idx, hap_t hap, bool matched = false) const;

    /** @brief Sets or unsets the alternate allele on one haplotype for a matched genotype. */
    void set_var_matched_gt_on_hap(int var_idx, hap_t hap, bool set = true,
            bool ignore_errors = false);

    /** @brief Classifies variant as SNP, INDEL, or SV based on reference length and g.sv_threshold. */
    sizeclass_t get_vartype(int vi);

    /** @brief Records a variant's allele count error type from its original and matched genotypes. */
    ac_errtype_t set_allele_errtype(int vi, bool query);

    /** @brief Returns true if haplotypes should be swapped when reporting matched_gt data relative to orig_gt. */
    bool matched_gt_is_swapped(int vi) const;

    // originally parsed data (size n)
    std::string ctg;                ///< Contig name (chromosome identifier)
    std::vector<int> poss;          ///< variant start positions (0-based)
    std::vector<int> rlens;         ///< reference lengths
    std::vector<edittype_t> types;     ///< variant type: NONE, SUB, INS, DEL, CPX
    std::vector<bedloc_t> locs;      ///< BED location: INSIDE, OUTSIDE, BORDER
    std::vector<std::string> refs;  ///< variant reference allele
    std::vector<std::string> alts;  ///< variant alternate allele (always one)
    std::vector<gt_t> orig_gts;  ///< simple genotype (0|1, 1|0, or 1|1)
    std::vector<float> gt_quals;    ///< genotype quality (capped above at --max-qual)
    std::vector<float> var_quals;   ///< variant quality (capped above at --max-qual)
    std::vector<int> phase_sets;    ///< integer representing variant phase set (0 = missing)
    std::vector<int> rec_idxs;      ///< source VCF record ordinal (0-based, -1 = unknown)
    std::vector<int> alt_idxs;      ///< original ALT ordinal (1-based, -1 = unknown)
    std::vector<ploidy_t> ploidies; ///< variant ploidy, from std::abs(ngt)
    std::vector<int> superclusters; ///< initially -1, set during superclustering
    int n = 0;                      ///< Total number of variants

    // set during clustering (size nc+1)
    std::vector<int> clusters;      ///< variant indices of clusters in this struct's vectors
    std::vector<int> left_reaches;  ///< cluster leftmost reach
    std::vector<int> right_reaches; ///< cluster rightmost reach
    int nc = 0;                     ///< Total number of clusters (size of clusters vector is nc+1 with sentinel)

    // set during prec_recall_aln() (size (2, n), additional axis for haplotype)
    std::vector<gt_t> matched_gts;  ///< the other callset's genotype (0|1, 1|0, or 1|1) recovered by alignment
    EnumArray<hap_t, std::vector<errtype_t>, HAP_SLOTS> errtypes; ///< error type: TP, FP, FN
    EnumArray<hap_t, std::vector<int>, HAP_SLOTS> sync_group;   ///< group of variants that participate in credit
    EnumArray<hap_t, std::vector<float>, HAP_SLOTS> callq;      ///< min call quality in sync group (for truth, of associated call)
    EnumArray<hap_t, std::vector<int>, HAP_SLOTS> ref_ed;       ///< reference edit distance in sync group
    EnumArray<hap_t, std::vector<int>, HAP_SLOTS> query_ed;     ///< query edit distance in sync group
    EnumArray<hap_t, std::vector<float>, HAP_SLOTS> credit;     ///< percentage reduction in edit dist (ref->query)

    // set during phase() (size n)
    std::vector<phase_t> phases;     ///< variant keep/swap/unknown, from alignment (matched_gt relative to orig_gt)
    std::vector<phase_t> pb_phases;  ///< phaseblock keep/swap, from phasing algorithm
    std::vector<ac_errtype_t> ac_errtype; ///< allele count error type, truth count then query count on both callsets (e.g. 0|1 -> 1|1)

    // shared with every other container of this callset, indexed by rec_idxs
    std::shared_ptr<srcRecords> src_recs; ///< retained source records (nullptr = none retained)
};

/**
 * @class variantData
 * @brief Store variant data across all haplotypes and contigs for a single callset.
 */
class variantData {
public:
    /** @brief Constructs an empty variant data container defaulting to QUERY callset. */
    variantData();

    // data
    std::shared_ptr<fastaData> ref;  ///< Pointer to reference FASTA data
    callset_t callset;               ///< Callset type: QUERY or TRUTH
    std::string filename;            ///< Source VCF filename

    std::string sample;               ///< Sample name from VCF header
    std::vector<std::string> contigs; ///< List of all contig names
    std::vector<int> lengths;         ///< List of all contig lengths
    std::vector<                      ///< Ploidies observed on each contig, parallel to contigs
        std::set<int> > observed_ploidies;
    ///< Per-haplotype, per-contig variant containers: variants[hap][ctg]
    EnumArray<hap_t,
        std::unordered_map<std::string, std::shared_ptr<ctgVariants> >, HAP_SLOTS> variants;
    std::unordered_map<                ///< Per-contig containers of retained, unevaluated variants
        std::string,
        std::shared_ptr<ctgSideline> > sidelined;
    std::shared_ptr<srcRecords>       ///< Source records retained from this callset's VCF
        src_recs;
};

/** @brief Classifies a record's raw GT array into its parse-time genotype shape. */
gtparse_t classify_gt(const int32_t * gt, int ngt);

/** @brief Builds the summary VCF header, declaring every FORMAT field and the TRUTH/QUERY samples. */
bcf_hdr_t* summary_vcf_header(const std::vector<std::string> & contigs,
        const std::vector<int> & lengths,
        const EnumArray<callset_t, std::shared_ptr<srcRecords>, CALLSET_SLOTS> & src_recs = {});

/** @brief Returns the FORMAT values of a sample that made no call at a locus. */
sample_fields empty_sample_fields(int sc_idx, int phase_block);

/** @brief Sets every FORMAT field of one record from the two samples' values. */
void set_record_samples(const bcf_hdr_t* hdr, bcf1_t* rec,
        const sample_fields & truth, const sample_fields & query);

/** @brief Parses variants from a VCF file into a variantData container, with filtering and validation. */
void parse_variants(const std::string & vcf_fn,
        std::shared_ptr<variantData> variant_data,
        std::shared_ptr<fastaData> reference, callset_t callset);

#endif
