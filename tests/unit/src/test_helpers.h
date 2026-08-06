/**
 * @file test_helpers.h
 * @brief Shared unit-test scaffolding: global-state fixture, temporary files, in-memory builders.
 */
#ifndef _TEST_HELPERS_H_
#define _TEST_HELPERS_H_

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../../../src/bed.h"
#include "../../../src/cluster.h"
#include "../../../src/defs.h"
#include "../../../src/dist.h"
#include "../../../src/fasta.h"
#include "../../../src/globals.h"
#include "../../../src/phase.h"
#include "../../../src/variant.h"

/* Global state fixture ***************************************************************************/

/**
 * @class GlobalsGuard
 * @brief Saves the global `g` on construction and restores it on destruction, silencing logging.
 *
 * Declare one at the top of every test that reads or writes `g`. Verbosity is set to 0 so that
 * INFO output does not pollute the test log. Restoration is a whole-struct assignment, so every
 * field is covered and a newly added `Globals` field needs no change here.
 */
class GlobalsGuard {
public:
    /** @brief Saves the current global configuration and sets verbosity to 0. */
    GlobalsGuard();

    /** @brief Restores the saved global configuration. */
    ~GlobalsGuard();

    GlobalsGuard(const GlobalsGuard &) = delete;
    GlobalsGuard & operator=(const GlobalsGuard &) = delete;

private:
    Globals saved; ///< Copy of the global configuration taken at construction
};

/* Temporary files ********************************************************************************/

/**
 * @class TempDir
 * @brief Creates a uniquely-named temporary directory and removes it recursively on destruction.
 */
class TempDir {
public:
    /** @brief Creates a uniquely-named directory beneath the system temporary directory. */
    TempDir(const std::string & prefix = "vcfdist_test");

    /** @brief Recursively removes the directory and everything beneath it. */
    ~TempDir();

    TempDir(const TempDir &) = delete;
    TempDir & operator=(const TempDir &) = delete;

    /** @brief Returns the directory path, without a trailing slash. */
    std::string path() const;

    /** @brief Returns the path of a named entry within the directory. */
    std::string path(const std::string & name) const;

private:
    std::string dir; ///< Absolute path of the created directory
};

/**
 * @struct vcf_opts
 * @brief Header lines and sample name used to write a temporary single-sample VCF.
 *
 * Each group of header lines is replaceable on its own, so a test can drop the contig `length`
 * key or omit a FILTER declaration without restating the rest of the header.
 */
struct vcf_opts {
    std::string filename = "test.vcf"; ///< Basename of the VCF written within the TempDir
    std::string sample = "TRUTH";      ///< Sample name in the #CHROM line
    std::vector<std::string> meta = {"##fileformat=VCFv4.2"}; ///< Leading meta-information lines
    std::vector<std::string> contigs = ///< ##contig lines
        {"##contig=<ID=chr1,length=12>"};
    std::vector<std::string> filters = ///< ##FILTER lines
        {"##FILTER=<ID=PASS,Description=\"All filters passed\">"};
    std::vector<std::string> infos;    ///< ##INFO lines, none by default
    std::vector<std::string> formats = ///< ##FORMAT lines
        {"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
         "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
         "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">"};
};

/**
 * @brief Writes a minimal single-sample VCF into a temporary directory and returns its path.
 *
 * The directory is passed in rather than created here so that the caller controls its lifetime;
 * the file is deleted along with the directory. Records are written verbatim, one line each.
 */
std::string write_tmp_vcf(const TempDir & dir, const std::vector<std::string> & records,
        const vcf_opts & opts = vcf_opts());

/**
 * @brief Writes a BED file into a temporary directory and returns its path.
 *
 * Lines are written verbatim, so a test can supply extra columns or malformed records that the
 * in-memory make_bed() builders cannot express.
 */
std::string write_tmp_bed(const TempDir & dir, const std::vector<std::string> & lines,
        const std::string & name = "test.bed");

/** @brief Returns the path of a checked-in fixture under tests/unit/data/. */
std::string data_path(const std::string & name);

/** @brief Reads an entire file into a string, yielding "" if it cannot be opened. */
std::string read_text(const std::string & fn);

/**
 * @brief Builds header options for one sample, declaring each named contig with one length.
 *
 * Subsumes the boilerplate of writing "##contig=<ID=...,length=...>" by hand. A test that needs a
 * malformed contig line, an extra FILTER, or a FORMAT declaration dropped assigns over the
 * corresponding vcf_opts field afterwards.
 */
vcf_opts make_vcf_opts(callset_t callset = QUERY,
        const std::vector<std::string> & contigs = {"chr1"}, int length = 1000);

/* Record lines ***********************************************************************************/

/**
 * @struct vcf_record
 * @brief Columns of one VCF data line, so a caller can vary only the column it targets.
 *
 * Every column is a string except POS, so a test can write "." for an unreported QUAL or a
 * comma-separated ALT list without a separate builder for each shape.
 */
struct vcf_record {
    std::string ctg = "chr1";      ///< CHROM column
    int pos = 100;                 ///< POS column, 1-based
    std::string id = ".";          ///< ID column
    std::string ref = "A";         ///< REF allele
    std::string alt = "G";         ///< ALT allele, comma-separated when multiallelic
    std::string qual = "50";       ///< QUAL column, "." when no quality is reported
    std::string filter = "PASS";   ///< FILTER column, "." when no filters are applied
    std::string info = ".";        ///< INFO column
    std::string format = "GT:PS";  ///< FORMAT column
    std::string sample = "1|0:1";  ///< Sample column, matching format field for field
};

/** @brief Joins the columns of one VCF data line with tabs, without a trailing newline. */
std::string vcf_line(const vcf_record & rec);

/** @brief Builds a data line with the given genotype, carrying QUAL 50, PASS, and phase set 1. */
std::string record(int pos, const std::string & ref, const std::string & alt,
        const std::string & gt, const std::string & ctg = "chr1");

/** @brief Builds a data line carrying the given FORMAT keys and sample values. */
std::string fmt_record(int pos, const std::string & ref, const std::string & alt,
        const std::string & format, const std::string & sample,
        const std::string & ctg = "chr1");

/** @brief Builds a phased 1|0 SNP with the given QUAL and FILTER columns. */
std::string qual_filter_record(int pos, const std::string & qual, const std::string & filter,
        const std::string & ctg = "chr1");

/* Parse capture **********************************************************************************/

/**
 * @class StderrToFile
 * @brief Redirects the C `stderr` stream to a file for the object's lifetime.
 *
 * WARN() and INFO() reach stderr through fprintf(), so capturing the summary means redirecting
 * the underlying file descriptor; swapping std::cerr's streambuf would not intercept it.
 * Restoring in the destructor keeps a failed assertion from leaving stderr pointing into the
 * temporary directory after TempDir has deleted it.
 */
class StderrToFile {
public:
    /** @brief Redirects stderr to fn, truncating any existing contents. */
    explicit StderrToFile(const std::string & fn);

    /** @brief Flushes the redirected output and restores the original stderr. */
    ~StderrToFile();

    StderrToFile(const StderrToFile &) = delete;
    StderrToFile & operator=(const StderrToFile &) = delete;

private:
    int saved_fd; ///< Duplicate of the original stderr descriptor
    int file_fd;  ///< Descriptor of the redirect target
};

/**
 * @struct ParseResult
 * @brief Parsed variants plus everything parse_variants() reported to stderr.
 */
struct ParseResult {
    std::shared_ptr<variantData> vars; ///< Variants that survived parse-time filtering
    std::string log;                   ///< All INFO/WARN output from parse_variants()
};

/** @brief Parses VCF records with parse_variants(), capturing its stderr. */
ParseResult parse_records(const TempDir & dir, const std::vector<std::string> & records);

/**
 * @brief Parses VCF records under a caller-supplied header, capturing stderr.
 *
 * Lets a test drop a FORMAT declaration, rename a contig, or add a second sample without also
 * restating the reference; parse_variants() only stores the reference pointer, so the default
 * nullptr suffices unless a test asserts on variantData::ref.
 */
ParseResult parse_records(const TempDir & dir, const std::vector<std::string> & records,
        const vcf_opts & opts, std::shared_ptr<fastaData> ref = nullptr);

/**
 * @brief Parses records under the given header without redirecting stderr.
 *
 * Death tests match the ERROR message on stderr, so the redirect parse_records() installs would
 * hide it; the parsed output is discarded, since a parse that errors never returns.
 */
void parse_unredirected(const TempDir & dir, const std::vector<std::string> & records,
        const vcf_opts & opts, callset_t callset = QUERY);

/** @brief Reports whether a captured log contains a substring. */
bool logged(const std::string & log, const std::string & text);

/** @brief Reports whether the captured log contains a substring. */
bool logged(const ParseResult & r, const std::string & text);

/**
 * @brief Returns the variants that survived parsing on one haplotype of a contig.
 *
 * Yields nullptr for a contig the parse never reached, rather than inserting an empty entry the
 * way operator[] would; a caller that dereferences the result should assert on it first.
 */
std::shared_ptr<ctgVariants> hap_vars(const ParseResult & r, hap_t hap,
        const std::string & ctg = "chr1");

/** @brief Counts variants that survived parsing on one haplotype of a contig. */
int kept_on_hap(const ParseResult & r, hap_t hap, const std::string & ctg = "chr1");

/** @brief Counts variants that survived parsing across both haplotypes of a contig. */
int total_kept(const ParseResult & r, const std::string & ctg = "chr1");

/** @brief Reports whether any variant survived parsing at a VCF position on a contig. */
bool kept_pos(const ParseResult & r, int pos, const std::string & ctg = "chr1");

/**
 * @brief Counts the variants that survived parsing at a VCF position on a contig.
 *
 * Positions are the 1-based VCF coordinates the records were written with, so a test asserts on
 * the same number it passed to record(); an INS/DEL is stored at its first altered base, one past
 * the anchor base the record is positioned on.
 */
size_t count_pos(const ParseResult & r, int pos, const std::string & ctg = "chr1");

/** @brief Returns the genotype-histogram line parse_variants() prints for a genotype and count. */
std::string gt_hist_line(gtparse_t gt, int count);

/** @brief Returns the variant-type line parse_variants() prints for a type and count. */
std::string type_hist_line(edittype_t type, int count);

/* In-memory builders *****************************************************************************/

/** @brief Builds a single-contig reference, bypassing the FASTA-parsing constructor. */
std::shared_ptr<fastaData> make_fasta(const std::string & ctg, const std::string & seq,
        bool uppercase = true);

/** @brief Builds a multi-contig reference from (contig, sequence) pairs, in order. */
std::shared_ptr<fastaData> make_fasta(
        const std::vector< std::pair<std::string, std::string> > & seqs, bool uppercase = true);

/** @brief Builds a bedData holding the given [start, stop) regions on one contig, in order. */
bedData make_bed(const std::string & ctg, const std::vector< std::pair<int, int> > & regions);

/**
 * @brief Builds a multi-contig bedData from (contig, regions) pairs, in order.
 *
 * Regions are appended exactly as given, so an invalid layout can be built and handed to check().
 */
bedData make_bed(const std::vector< std::pair<std::string,
        std::vector< std::pair<int, int> > > > & regions);

/**
 * @brief Builds a variantData over the given contigs, each holding an empty ctgVariants per hap.
 *
 * Stands in for a parsed VCF in tests that care about the contig, length, and observed-ploidy
 * fields rather than about variants; `filename` and `sample` follow from the callset.
 */
std::shared_ptr<variantData> make_variantData(callset_t callset,
        const std::vector<std::string> & contigs, const std::vector<int> & lengths,
        const std::vector< std::set<int> > & observed_ploidies);

/**
 * @struct var_desc
 * @brief Describes one variant to append to a ctgVariants container.
 *
 * New fields are appended, never inserted, so that positional brace-initialization in existing
 * tests keeps binding to the same members.
 */
struct var_desc {
    int pos = 0;               ///< 0-based reference start position
    int rlen = 0;              ///< Reference allele length
    edittype_t type = TYPE_SUB; ///< Variant type (TYPE_SUB, TYPE_INS, TYPE_DEL, TYPE_CPX)
    std::string ref;           ///< Reference allele sequence
    std::string alt;           ///< Alternate allele sequence
    gt_t gt = GT_REF_ALT;     ///< Original genotype (GT_*)
    float qual = 60;           ///< Sets both var_qual and gt_qual (each clamped to g.max_qual)
    int phase_set = 0;         ///< Phase set identifier (0 = missing)
    int supercluster = -1;     ///< Supercluster index (-1 = not yet assigned)
    bedloc_t loc = BED_INSIDE; ///< BED location (BED_INSIDE, BED_OUTSIDE, BED_BORDER, BED_OFFCTG)
    int rec_idx = -1;          ///< Source VCF record ordinal, 0-based (-1 = unknown)
    int alt_idx = -1;          ///< Original ALT ordinal, 1-based (-1 = unknown)
    ploidy_t ploidy = PLOIDY_DIPLOID; ///< Variant ploidy
};

/** @brief Builds a ctgVariants container holding the described variants, in the given order. */
std::shared_ptr<ctgVariants> make_ctgVariants(const std::string & ctg,
        const std::vector<var_desc> & vars);

/**
 * @brief Builds a one-variant container with the given original and matched genotypes.
 *
 * The variant is an A>C substitution, since the genotype rather than the allele is what a caller
 * of this builder is varying.
 */
std::shared_ptr<ctgVariants> make_gt_var(gt_t orig_gt, gt_t matched_gt,
        const std::string & ctg = "chr1", int pos = 100);

/** @brief Builds a one-variant container of the given type with the given allele sequences. */
std::shared_ptr<ctgVariants> make_typed_var(edittype_t type, const std::string & ref,
        const std::string & alt, const std::string & ctg = "chr1", int pos = 100);

/** @brief Sets all six per-haplotype evaluation lanes for one variant. */
void set_hap_data(std::shared_ptr<ctgVariants> vars, hap_t hap, int idx, errtype_t errtype,
        int sync_group, float callq, int ref_ed, int query_ed, float credit);

/** @brief Sets cluster boundaries and reaches; nc defaults to clusters.size()-1. */
void set_clusters(std::shared_ptr<ctgVariants> vars, const std::vector<int> & clusters,
        const std::vector<int> & left_reaches, const std::vector<int> & right_reaches,
        int nc = -1);

/** @brief Builds a ctgSuperclusters holding the given query and truth variant containers. */
std::shared_ptr<ctgSuperclusters> make_ctgSuperclusters(std::shared_ptr<ctgVariants> qvars,
        std::shared_ptr<ctgVariants> tvars);

/** @brief Builds a superclusterData over the given contigs, bypassing clustering. */
std::shared_ptr<superclusterData> make_superclusterData(
        const std::vector<std::string> & contigs, const std::vector<int> & lengths,
        const std::vector< std::shared_ptr<ctgSuperclusters> > & superclusters,
        std::shared_ptr<fastaData> ref = nullptr);

/**
 * @brief Builds a phaseblockData over the given contigs, bypassing the phasing pipeline.
 *
 * The real constructor runs fix_phase_set_tags(), phase(), and fix_allele_counts(), each of which
 * would overwrite the matched genotypes, error types, and allele count error types a test set
 * by hand. It is therefore invoked over empty per-contig containers, and the caller's populated
 * ones are substituted afterwards.
 */
std::unique_ptr<phaseblockData> make_phaseblockData(
        const std::vector<std::string> & contigs, const std::vector<int> & lengths,
        const std::vector< std::shared_ptr<ctgSuperclusters> > & superclusters,
        std::shared_ptr<fastaData> ref = nullptr);

/** @brief Builds an alignment graph for one supercluster and truth haplotype. */
std::shared_ptr<Graph> make_graph(std::shared_ptr<ctgSuperclusters> sc,
        std::shared_ptr<fastaData> ref, const std::string & ctg, hap_t truth_hap, int sc_idx = 0);

/** @brief Allocates the offsets buffer that wf_swg_max_reach requires from its caller. */
std::vector<int> alloc_reach_offs(int qlen, int tlen, int x, int o, int e);

#endif
