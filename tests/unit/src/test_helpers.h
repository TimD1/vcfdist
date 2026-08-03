/**
 * @file test_helpers.h
 * @brief Shared unit-test scaffolding: global-state fixture, temporary files, in-memory builders.
 */
#ifndef _TEST_HELPERS_H_
#define _TEST_HELPERS_H_

#include <memory>
#include <string>
#include <vector>

#include "../../../src/cluster.h"
#include "../../../src/defs.h"
#include "../../../src/dist.h"
#include "../../../src/fasta.h"
#include "../../../src/globals.h"
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

/** @brief Returns the path of a checked-in fixture under tests/unit/data/. */
std::string data_path(const std::string & name);

/* In-memory builders *****************************************************************************/

/** @brief Builds a single-contig reference, bypassing the FASTA-parsing constructor. */
std::shared_ptr<fastaData> make_fasta(const std::string & ctg, const std::string & seq,
        bool uppercase = true);

/** @brief Builds a multi-contig reference from (contig, sequence) pairs, in order. */
std::shared_ptr<fastaData> make_fasta(
        const std::vector< std::pair<std::string, std::string> > & seqs, bool uppercase = true);

/**
 * @struct var_desc
 * @brief Describes one variant to append to a ctgVariants container.
 */
struct var_desc {
    int pos = 0;               ///< 0-based reference start position
    int rlen = 0;              ///< Reference allele length
    uint8_t type = TYPE_SUB;   ///< Variant type (TYPE_SUB, TYPE_INS, TYPE_DEL, TYPE_CPX)
    std::string ref;           ///< Reference allele sequence
    std::string alt;           ///< Alternate allele sequence
    uint8_t gt = GT_REF_ALT1;  ///< Original genotype (GT_*)
    float qual = 60;           ///< Variant and genotype quality (add_var clamps to g.max_qual)
    int phase_set = 0;         ///< Phase set identifier (0 = missing)
    int supercluster = -1;     ///< Supercluster index (-1 = not yet assigned)
    uint8_t loc = BED_INSIDE;  ///< BED location (BED_INSIDE, BED_OUTSIDE, BED_BORDER, BED_OFFCTG)
};

/** @brief Builds a ctgVariants container holding the described variants, in the given order. */
std::shared_ptr<ctgVariants> make_ctgVariants(const std::string & ctg,
        const std::vector<var_desc> & vars);

/** @brief Sets all six per-haplotype evaluation lanes for one variant. */
void set_hap_data(std::shared_ptr<ctgVariants> vars, int hap, int idx, uint8_t errtype,
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
        const std::vector<int> & ploidy,
        const std::vector< std::shared_ptr<ctgSuperclusters> > & superclusters,
        std::shared_ptr<fastaData> ref = nullptr);

/** @brief Builds an alignment graph for one supercluster and truth haplotype. */
std::shared_ptr<Graph> make_graph(std::shared_ptr<ctgSuperclusters> sc,
        std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hap, int sc_idx = 0);

/** @brief Allocates the offsets buffer that wf_swg_max_reach requires from its caller. */
std::vector<int> alloc_reach_offs(int qlen, int tlen, int x, int o, int e);

#endif
