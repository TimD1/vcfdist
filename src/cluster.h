#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "vcf.h"

#define PHASE_ORIG 0
#define PHASE_SWAP 1
#define PHASE_NONE 2

#define PHASE_PTR_KEEP 0
#define PHASE_PTR_SWAP 1


class ctgClusters {
public:

    // constructor
    ctgClusters() {;}

    // helper functions
    void add(
           variantCalls* cal1_vars, int cal1_beg_idx, int cal1_end_idx, 
           variantCalls* cal2_vars, int cal2_beg_idx, int cal2_end_idx, 
           variantCalls* hap1_vars, int hap1_beg_idx, int hap1_end_idx, 
           variantCalls* hap2_vars, int hap2_beg_idx, int hap2_end_idx, 
           int phase, int orig_phase_dist, int swap_phase_dist);

    variantCalls* cal1_vars = nullptr; 
    variantCalls* cal2_vars = nullptr; 
    variantCalls* hap1_vars = nullptr;
    variantCalls* hap2_vars = nullptr;
    std::vector<int> cal1_beg_idx, cal1_end_idx, cal2_beg_idx, cal2_end_idx,
        hap1_beg_idx, hap1_end_idx, hap2_beg_idx, hap2_end_idx;

    std::vector<int> phase;
    std::vector<int> orig_phase_dist, swap_phase_dist;
    int n = 0;
};

class clusterData {
public:
    clusterData(std::vector<std::string> ctgs) {
        for (auto ctg : ctgs) {
            this->contigs.push_back(ctg);
            clusters[ctg] = ctgClusters();
        }
    }
    std::vector<std::string> contigs;
    std::unordered_map<std::string, ctgClusters> clusters;
};

#endif
