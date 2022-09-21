#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <string>
#include <unordered_map>
#include <vector>

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
    void set_variants(
           std::shared_ptr<ctgVariants> calls1_vars,
           std::shared_ptr<ctgVariants> calls2_vars,
           std::shared_ptr<ctgVariants> truth1_vars,
           std::shared_ptr<ctgVariants> truth2_vars);
    void add_supercluster(
           int calls1_beg_idx, int calls1_end_idx, 
           int calls2_beg_idx, int calls2_end_idx, 
           int truth1_beg_idx, int truth1_end_idx, 
           int truth2_beg_idx, int truth2_end_idx, 
           int beg, int end,
           int phase, int orig_phase_dist, int swap_phase_dist);

    // pointers to variant/cluster data
    std::shared_ptr<ctgVariants> calls1_vars = nullptr; 
    std::shared_ptr<ctgVariants> calls2_vars = nullptr; 
    std::shared_ptr<ctgVariants> truth1_vars = nullptr;
    std::shared_ptr<ctgVariants> truth2_vars = nullptr;

    // cluster indices of superclusters
    std::vector<int> calls1_beg_idx, calls1_end_idx, calls2_beg_idx, calls2_end_idx,
        truth1_beg_idx, truth1_end_idx, truth2_beg_idx, truth2_end_idx;
    std::vector<int> begs, ends; // ref positions

    std::vector<int> phase;
    std::vector<int> orig_phase_dist, swap_phase_dist;
    int n = 0;
};

class clusterData {
public:
    clusterData(std::vector<std::string> ctgs) {
        for (auto ctg : ctgs) {
            this->contigs.push_back(ctg);
            ctg_superclusters[ctg] = std::shared_ptr<ctgClusters>(new ctgClusters());
        }
    }
    std::vector<std::string> contigs;
    std::unordered_map<std::string, std::shared_ptr<ctgClusters> > ctg_superclusters;
};

void cluster(std::unique_ptr<variantData> & vcf);

#endif
