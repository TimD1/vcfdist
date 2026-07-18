/**
 * @file globals.h
 * @brief Global configuration object and program-wide string table declarations.
 */
#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "bed.h"
#include "defs.h"
#include "timer.h"

/**
 * @class Globals
 * @brief Global configuration and program parameters.
 *
 * Single instance `g` is accessible throughout the program and populated by `parse_args()`.
 */
class Globals {
public:
    // constructors
    Globals() {;}

    // input files
    std::string ref_fasta_fn; ///< Reference FASTA filename
    FILE* ref_fasta_fp;       ///< Open file pointer to reference FASTA
    std::string query_vcf_fn; ///< Query VCF filename
    std::string truth_vcf_fn; ///< Truth VCF filename
    std::string bed_fn;       ///< BED file filename (empty if not provided)
    bedData bed;              ///< Parsed BED region data
    bool bed_exists = false;  ///< True if a BED file was provided
    bool write = true;        ///< If false, skip writing output files

    // variant filtering
    std::vector<std::string> filters; ///< FILTER tag values to include (empty = all)
    std::vector<int> filter_ids;      ///< htslib integer IDs corresponding to filters
    int min_qual = 0;                 ///< Minimum variant quality score (inclusive)
    int max_qual = 60;                ///< Maximum variant quality score (variants above are clamped)
    int max_size = 1000;              ///< Maximum variant size in bp; larger variants are ignored

    // clustering
    int max_supercluster_size = 15000; ///< Maximum allowed supercluster size in bp before splitting
    int cluster_min_gap = 50;          ///< Minimum gap between clusters (for gap/size methods)
    int reach_min_gap = 10;            ///< Minimum gap for reach-based cluster merging
    int max_cluster_itrs = 1;          ///< Maximum iterations for expanding/merging clusters
    int sub = 5;                       ///< Smith-Waterman substitution penalty
    int open = 6;                      ///< Smith-Waterman gap-open penalty
    int extend = 2;                    ///< Smith-Waterman gap-extend penalty

    // precision-recall
    double credit_threshold = 0.98; ///< Minimum partial credit to count a variant as TP
    int max_dist = 100;             ///< Maximum alignment edit distance per supercluster
    int max_retries = 0;            ///< Maximum retries for alignment after removing large variants

    // memory params
    int max_threads = 64;          ///< Maximum number of threads to use
    double max_ram = 64;           ///< Approximate maximum RAM in GB for alignment
    int thread_nsteps;             ///< Number of thread/RAM scheduling steps
    std::vector<int> thread_steps; ///< Thread counts at each scheduling step
    std::vector<float> ram_steps;  ///< RAM-per-thread values at each scheduling step

    // high-level options
    int verbosity = 1;      ///< Verbosity level: 0 = succinct, 1 = default, 2 = verbose
    int sv_threshold = 50;  ///< Variants >= this size (bp) are classified as SVs
    std::string out_prefix; ///< Output file path prefix
    std::string cmd;        ///< Full command-line string as provided

    // member functions
    /** @brief Parses command-line arguments and initializes global configuration. */
    void parse_args(int argc, char ** argv);

    /** @brief Prints program version to stdout. */
    void print_version() const;

    /** @brief Prints usage information and all command-line options to stdout. */
    void print_usage() const;

    /** @brief Prints publication citation in MLA and BibTeX formats. */
    void print_citation() const;

    /** @brief Initializes one named timer object per pipeline stage. */
    void init_timers(const std::vector<std::string> & timer_strs);

    // program data
    const std::string VERSION = "3.0.0-b0"; ///< Program version string
    const std::string PROGRAM = "vcfdist";  ///< Program name string
    std::vector<timer> timers;              ///< Per-stage pipeline timers
};

/** @brief Extracts the parent directory path from a file path string. */
std::string parent_path(const std::string & out_prefix);

/** @brief Creates a directory and all necessary parent directories. */
void create_directory(const std::string & dir);

extern Globals g; ///< Global program configuration instance

// defined in main.cpp
extern std::vector<std::string> callset_strs; ///< String representations of QUERY/TRUTH callset indices
extern std::vector<std::string> error_strs;   ///< String representations of ERRTYPE_* constants
extern std::vector<std::string> gt_strs;      ///< String representations of GT_* genotype constants
extern std::vector<std::string> ac_strs;      ///< String representations of AC_ERR_* allele count error types
extern std::vector<std::string> phase_strs;   ///< String representations of PHASE_* constants
extern std::vector<std::string> region_strs;  ///< String representations of BED_* location constants
extern std::vector<std::string> switch_strs;  ///< String representations of SWITCHTYPE_* constants
extern std::vector<std::string> timer_strs;   ///< String names for pipeline stage timers (TIME_* order)
extern std::vector<std::string> type_strs;    ///< String representations of TYPE_* variant type constants
extern std::vector<std::string> type_strs2;   ///< Alternate string representations of TYPE_* constants (ALL/SNP/INS/DEL/INDEL)
extern std::vector<std::string> vartype_strs; ///< String representations of VARTYPE_* size-class constants

#endif
