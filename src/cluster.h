#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "variant.h"

#define PHASE_ORIG 0
#define PHASE_SWAP 1
#define PHASE_NONE 2

#define PHASE_PTR_KEEP 0
#define PHASE_PTR_SWAP 1


class ctgClusters {
public:

    // constructor
    ctgClusters() {;}

    // add supercluster info
    void add_supercluster(
           int query1_beg_idx, int query1_end_idx, 
           int query2_beg_idx, int query2_end_idx, 
           int truth1_beg_idx, int truth1_end_idx, 
           int truth2_beg_idx, int truth2_end_idx, 
           int beg, int end);
    void add_phasing(
           int phase, 
           int orig_phase_dist, 
           int swap_phase_dist);

    // pointers to variant data
    std::shared_ptr<ctgVariants> query1_vars = nullptr; 
    std::shared_ptr<ctgVariants> query2_vars = nullptr; 
    std::shared_ptr<ctgVariants> truth1_vars = nullptr;
    std::shared_ptr<ctgVariants> truth2_vars = nullptr;

    // cluster indices of superclusters
    std::vector<int> query1_beg_idx, query1_end_idx, query2_beg_idx, query2_end_idx,
        truth1_beg_idx, truth1_end_idx, truth2_beg_idx, truth2_end_idx;
    std::vector<int> begs, ends; // ref positions
    int n = 0; // number of superclusters

    // phasing information per supercluster
    std::vector<int> phase;
    std::vector<int> orig_phase_dist, swap_phase_dist;
};

class clusterData {
public:
    clusterData(
            std::unique_ptr<variantData> & query_ptr,
            std::unique_ptr<variantData> & truth_ptr,
            std::shared_ptr<fastaData> ref_ptr);

    void gap_supercluster();

    // data
    std::vector<std::string> contigs;
    std::unordered_map<std::string, std::shared_ptr<ctgClusters> > ctg_superclusters;
    std::shared_ptr<fastaData> ref;
};


// for single haplotype clustering (one VCF)
void cluster(std::unique_ptr<variantData> & vcf);
void sw_cluster(std::unique_ptr<variantData> & vcf, int sub, int open, int extend);

#endif
