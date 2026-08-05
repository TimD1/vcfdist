/**
 * @file variant.h
 * @brief Per-contig and per-callset variant containers with VCF parsing and output utilities.
 */
#ifndef _VARIANT_H_
#define _VARIANT_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#include "htslib/vcf.h"

#include "fasta.h"
#include "defs.h"

/**
 * @class ctgVariants
 * @brief Store all variant information for a single contig and callset.
 */
class ctgVariants {
public:

    /** @brief Constructs a contig-specific variant container with empty data structures. */
    ctgVariants(const std::string & ctg);

    /** @brief Appends a variant with all fields explicitly specified. */
    void add_var(int pos, int rlen, uint8_t type, uint8_t loc,
        const std::string & ref, const std::string & alt, uint8_t orig_gt, float gt_qual, float var_qual,
        int phase_set, int rec_idx = -1, int alt_idx = -1, uint8_t ploidy = 0,
        int supercluster = -1, uint8_t calc_gt = GT_REF_REF,
        uint8_t hap1_errtype = ERRTYPE_UN, uint8_t hap2_errtype = ERRTYPE_UN,
        int hap1_sync_group = 0, int hap2_sync_group = 0, float hap1_callq = 0, float hap2_callq = 0,
        int hap1_ref_ed = 0, int hap2_ref_ed = 0, int hap1_query_ed = 0, int hap2_query_ed = 0,
        float hap1_credit = 0, float hap2_credit = 0);

    /** @brief Writes fixed VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) for one variant. */
    void print_var_info(FILE* out_fp, std::shared_ptr<fastaData> ref,
            const std::string & ctg, int idx);

    /** @brief Writes dot-separated empty sample fields for a variant with no call on this haplotype. */
    void print_var_empty(FILE* out_fp, int sc_idx, int phase_block, bool query = false);

    /** @brief Writes sample-specific FORMAT fields for one variant to output VCF. */
    void print_var_sample(FILE* out_fp, int vi, int hi, const std::string & gt, int sc_idx,
            int phase_block, bool phase_switch, bool phase_flip, bool query = false);

    /** @brief Returns true if a variant is present on the specified haplotype. */
    bool var_on_hap(int var_idx, int hap, bool calc = false) const;

    /** @brief Sets or unsets the alternate allele on one haplotype for a calculated genotype. */
    void set_var_calcgt_on_hap(int var_idx, int hap, bool set = true, bool ignore_errors = false);

    /** @brief Classifies variant as SNP, INDEL, or SV based on reference length and g.sv_threshold. */
    int get_vartype(int vi);

    /** @brief Records allele count error type by comparing original and calculated genotypes. */
    int set_allele_errtype(int vi);

    /** @brief Returns true if haplotypes should be swapped when reporting calc_gt data relative to orig_gt. */
    bool calcgt_is_swapped(int vi) const;

    // originally parsed data (size n)
    std::string ctg;                ///< Contig name (chromosome identifier)
    std::vector<int> poss;          ///< variant start positions (0-based)
    std::vector<int> rlens;         ///< reference lengths
    std::vector<uint8_t> types;     ///< variant type: NONE, SUB, INS, DEL, CPX
    std::vector<uint8_t> locs;      ///< BED location: INSIDE, OUTSIDE, BORDER
    std::vector<std::string> refs;  ///< variant reference allele
    std::vector<std::string> alts;  ///< variant alternate allele (always one)
    std::vector<uint8_t> orig_gts;  ///< simple genotype (0|1, 1|0, or 1|1)
    std::vector<float> gt_quals;    ///< genotype quality (capped above at --max-qual)
    std::vector<float> var_quals;   ///< variant quality (capped above at --max-qual)
    std::vector<int> phase_sets;    ///< integer representing variant phase set (0 = missing)
    std::vector<int> rec_idxs;      ///< source VCF record ordinal (0-based, -1 = unknown)
    std::vector<int> alt_idxs;      ///< original ALT ordinal (1-based, -1 = unknown)
    std::vector<uint8_t> ploidies;  ///< variant ploidy from std::abs(ngt) (0 = unknown)
    std::vector<int> superclusters; ///< initially -1, set during superclustering
    int n = 0;                      ///< Total number of variants

    // set during clustering (size nc+1)
    std::vector<int> clusters;      ///< variant indices of clusters in this struct's vectors
    std::vector<int> left_reaches;  ///< cluster leftmost reach
    std::vector<int> right_reaches; ///< cluster rightmost reach
    int nc = 0;                     ///< Total number of clusters (size of clusters vector is nc+1 with sentinel)

    // set during prec_recall_aln() (size (2, n), additional axis for haplotype)
    std::vector<uint8_t> calc_gts;                ///< calculated genotype (0|1, 1|0, or 1|1), only set for query
    std::vector< std::vector<uint8_t> > errtypes; ///< error type: TP, FP, FN
    std::vector< std::vector<int> > sync_group;   ///< group of variants that participate in credit
    std::vector< std::vector<float> > callq;      ///< min call quality in sync group (for truth, of associated call)
    std::vector< std::vector<int> > ref_ed;       ///< reference edit distance in sync group
    std::vector< std::vector<int> > query_ed;     ///< query edit distance in sync group
    std::vector< std::vector<float> > credit;     ///< percentage reduction in edit dist (ref->query)

    // set during phase() (size n)
    std::vector<int> phases;     ///< variant keep/swap/unknown, from alignment (calc_gt relative to orig_gt)
    std::vector<int> pb_phases;  ///< phaseblock keep/swap, from phasing algorithm
    std::vector<int> ac_errtype; ///< allele count error type (e.g. 0|1 -> 1|1)
};

/**
 * @class variantData
 * @brief Store variant data across all haplotypes and contigs for a single callset.
 */
class variantData {
public:
    /** @brief Constructs an empty variant data container defaulting to QUERY callset. */
    variantData();

    /** @brief Writes all parsed variants to a phased VCF file. */
    void write_vcf(std::string out_vcf_fn);

    /** @brief Writes a single variant record to a VCF file with GT and PS FORMAT fields. */
    void print_variant(FILE* out_fp, const std::string & ctg, int pos, int type,
        const std::string & ref, const std::string & alt, float qual, int phase_set,
        const std::string & gt);

    /** @brief Parses a CIGAR string and adds resulting variants to the container. */
    void add_variants(const std::vector<int> & cigar, int hap,
            int ref_pos, const std::string & ctg, const std::string & query,
            const std::string & ref, int qual, int phase_set);

    // data
    std::shared_ptr<fastaData> ref;  ///< Pointer to reference FASTA data
    int callset;                     ///< Callset type: QUERY (0) or TRUTH (1)
    std::string filename;            ///< Source VCF filename

    std::string sample;               ///< Sample name from VCF header
    std::vector<std::string> contigs; ///< List of all contig names
    std::vector<int> lengths;         ///< List of all contig lengths
    std::vector<int> ploidy;          ///< List of all contig ploidies
    std::vector<                      ///< Per-haplotype, per-contig variant containers: variants[hap][ctg]
        std::unordered_map<
            std::string,
            std::shared_ptr<ctgVariants> > > variants;
};

/** @brief Parses variants from a VCF file into a variantData container, with filtering and validation. */
void parse_variants(const std::string & vcf_fn,
        std::shared_ptr<variantData> variant_data,
        std::shared_ptr<fastaData> reference, int callset);

#endif
