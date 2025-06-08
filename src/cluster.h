#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "fasta.h"
#include "variant.h"
#include "defs.h"

class ctgSuperclusters {
public:
    ctgSuperclusters();

    // stores variant info for each callset
    // callset_vars[truth|query]
    std::vector< std::shared_ptr<ctgVariants> > callset_vars;

    int get_min_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end);
    int get_max_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end);
};

class superclusterData {
public:
    superclusterData(
            std::shared_ptr<variantData> query_ptr,
            std::shared_ptr<variantData> truth_ptr,
            std::shared_ptr<fastaData> ref_ptr);

    void add_callset_vars(int callset, std::vector< std::unordered_map< std::string,
            std::shared_ptr<ctgVariants> > > & vars);
    void supercluster(bool print = false);

    // data
    std::vector<std::string> contigs;
    std::vector<std::string> samples;
    std::vector<std::string> filenames;
    std::vector<int> lengths;
    std::vector<int> ploidy;
    std::unordered_map<std::string, 
        std::shared_ptr<ctgSuperclusters> > superclusters;
    std::shared_ptr<fastaData> ref;
};


// for single haplotype clustering (one VCF)
void simple_cluster(std::shared_ptr<variantData> vcf, int callset);
void wf_swg_cluster(variantData * vcf, int ctg_idx, int hap,
        int sub, int open, int extend);
std::vector< std::vector< std::vector<int> > > 
        sort_superclusters(std::shared_ptr<superclusterData>);
void duplicate_variant_info(std::shared_ptr<superclusterData>);


/******************************************************************************/


/* Supercluster splits are considered directly prior to each variant. The variant is at position
 * `start_pos` to `end_pos` on callset `callset_idx`.
 */
struct var_info {
    int callset_idx = 0;
    int start_pos = 0;
    int end_pos = 0;

    var_info(int _callset_idx, int _start_pos, int _end_pos) {
        this->callset_idx = _callset_idx;
        this->start_pos = _start_pos;
        this->end_pos = _end_pos;
    }
};

// helper functions for splitting large superclusters
std::vector< std::vector<int> > split_large_supercluster(
        std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        std::vector<int> & cluster_end_indices, bool print = false);
std::vector<int> get_supercluster_range(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices);
std::vector<int> get_supercluster_split_location(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices, bool print = false);
std::vector<int> split_cluster(
        std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & variant_split_indices,
        std::vector< std::vector<int> > & breakpoints, int breakpoint_idx, bool print = false);
var_info get_next_variant_info(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & var_curr_indices,
        const std::vector<int> & var_end_indices);

#endif
