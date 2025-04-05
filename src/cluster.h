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

    // constructor
    ctgSuperclusters();

    // add supercluster info
    void add_supercluster(
           std::vector<int> brks,
           int beg, int end);
    void set_phase(
           int sc_idx, int phase, 
           int orig_phase_dist, 
           int swap_phase_dist);

    // stores variant info for each contig
    // ctg_variants[truth/query][hap]
    std::vector< std::vector< std::shared_ptr<ctgVariants> > > ctg_variants;

    // cluster indices of superclusters (not variant indices!)
    // superclusters[truth/query][hap] = n+1
    std::vector< std::vector< std::vector<int> > > superclusters;

    // data (length n)
    int n = 0;                                         // number of superclusters
    std::vector<int> begs, ends;                       // ref beg/end pos across haps and truth/query
    std::vector<int> sc_phase;                         // keep/swap/unknown, from alignment
    std::vector<int> pb_phase;                         // keep/swap, from phasing algorithm
    std::vector<int> phase_sets;                       // input, from first variant in sc
    std::vector<int> orig_phase_dist, swap_phase_dist; // alignment distances
};

class superclusterData {
public:
    superclusterData(
            std::shared_ptr<variantData> query_ptr,
            std::shared_ptr<variantData> truth_ptr,
            std::shared_ptr<fastaData> ref_ptr);

    void supercluster(bool print = false);
    void transfer_phase_sets();

    // data
    std::vector<std::string> contigs;
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

/******************************************************************************/

/* Supercluster splits are considered directly prior to each variant. The variant is at position
 * `pos` on callset `hap_idx>>1` and hap `hap_idx&1`, with gap `gap` between the previous variant.
 */
struct var_info {
    int hap_idx = 0;
    int start_pos = 0;
    int end_pos = 0;

    var_info(int _hap_idx, int _start_pos, int _end_pos) {
        this->hap_idx = _hap_idx;
        this->start_pos = _start_pos;
        this->end_pos = _end_pos;
    }
};

// helper functions for splitting large superclusters
std::vector< std::vector<int> > split_large_supercluster(
        std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & cluster_start_indices,
        std::vector<int> & cluster_end_indices, bool print = false);
std::vector<int> get_supercluster_range(
        const std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices);
std::vector<int> get_supercluster_split_location(
        const std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices, bool print = false);
std::vector<int> split_cluster(
        std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & variant_split_indices,
        std::vector< std::vector<int> > & breakpoints, int breakpoint_idx, bool print = false);
var_info get_next_variant_info(
        const std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & var_curr_indices,
        const std::vector<int> & var_end_indices);

#endif
