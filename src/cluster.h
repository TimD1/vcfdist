#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "fasta.h"
#include "variant.h"
#include "defs.h"

/**
 * @class ctgSuperclusters
 * @brief Store truth and query variants from a given contig.
 */
class ctgSuperclusters {
public:
    /** @brief Constructs an empty contig supercluster container. */
    ctgSuperclusters();

    /** @brief ctgVariants info for each callset, indexed by TRUTH or QUERY */
    std::vector< std::shared_ptr<ctgVariants> > callset_vars;

    /** @brief Returns minimum reference position across a range of query and truth variants. */
    int get_min_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end);

    /** @brief Returns maximum reference position across a range of query and truth variants. */
    int get_max_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end);
};

/**************************************************************************************************/

/**
 * @class superclusterData
 * @brief Store superclustered variant data.
 */
class superclusterData {
public:
    /** @brief Constructs supercluster container from query and truth variant data and reference. */
    superclusterData(
            std::shared_ptr<variantData> query_ptr,
            std::shared_ptr<variantData> truth_ptr,
            std::shared_ptr<fastaData> ref_ptr);

    /** @brief Merges per-haplotype variant data for one callset into single per-contig containers. */
    void load_and_merge_callset_vars_across_haps(int callset,
            std::vector< std::unordered_map< std::string, std::shared_ptr<ctgVariants> > > & vars);

    /** @brief Groups variants into superclusters where truth and query variants may interact. */
    void supercluster(bool print = false);

    std::vector<std::string> samples;    ///< list containing QUERY and TRUTH VCF SAMPLE names
    std::vector<std::string> filenames;  ///< list containing QUERY and TRUTH VCF filenames
    std::vector<std::string> contigs;    ///< list of all contig names
    std::vector<int> lengths;            ///< list of all contig lengths
    std::vector<int> ploidy;             ///< list of all contig ploidies
    std::unordered_map<std::string,      ///< map from contig names to ctgSuperclusters
        std::shared_ptr<ctgSuperclusters> > superclusters;
    std::shared_ptr<fastaData> ref;      ///< pointer to reference fastaData
};

/**************************************************************************************************/

/**
 * @struct var_info
 * @brief Stores the interval of a variant on the reference.
 *
 * Supercluster splits are considered directly prior to each variant.
 */
struct var_info {
    int callset_idx = 0; ///< whether the variant is on the TRUTH or QUERY
    int start_pos = 0;   ///< the 0-based inclusive start position of the variant
    int end_pos = 0;     ///< the 0-based exclusive end position of the variant

    /** @brief Constructs a variant interval on the reference. */
    var_info(int _callset_idx, int _start_pos, int _end_pos) {
        this->callset_idx = _callset_idx;
        this->start_pos = _start_pos;
        this->end_pos = _end_pos;
    }
};

/**************************************************************************************************/

/** @brief Clusters variants on a single haplotype using gap-based methods. */
void simple_cluster(std::shared_ptr<variantData> vcf, int callset);

/** @brief Clusters variants using wavefront Smith-Waterman alignment. */
void wf_swg_cluster(variantData * vcf, int ctg_idx, int hap,
        int sub, int open, int extend);

/** @brief Returns superclusters sorted by size for multi-threaded scheduling. */
std::vector< std::vector< std::vector<int> > >
        sort_superclusters(std::shared_ptr<superclusterData>);

/***************************************************************************************************/

/** @brief Splits an oversized supercluster into smaller pieces at optimal breakpoints. */
std::vector< std::vector<int> > split_large_supercluster(
        std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        std::vector<int> & cluster_end_indices, bool print = false);

/** @brief Returns [beg_pos, end_pos] genomic range covered by a set of clusters. */
std::vector<int> get_supercluster_range(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices);

/** @brief Identifies optimal variant indices at which to split a supercluster. */
std::vector<int> get_supercluster_split_location(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices, bool print = false);

/** @brief Splits cluster boundaries at a given variant index and returns new cluster indices. */
std::vector<int> split_cluster(
        std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & variant_split_indices,
        std::vector< std::vector<int> > & breakpoints, int breakpoint_idx, bool print = false);

/** @brief Returns interval of the next unprocessed variant across both callsets. */
var_info get_next_variant_info(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & var_curr_indices,
        const std::vector<int> & var_end_indices);

#endif
