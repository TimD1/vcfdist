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

    // add supercluster info
    void set_variants(
           std::shared_ptr<variantCalls> cal1_vars,
           std::shared_ptr<variantCalls> cal2_vars,
           std::shared_ptr<variantCalls> hap1_vars,
           std::shared_ptr<variantCalls> hap2_vars);
    void add_supercluster(
           int cal1_beg_idx, int cal1_end_idx, 
           int cal2_beg_idx, int cal2_end_idx, 
           int hap1_beg_idx, int hap1_end_idx, 
           int hap2_beg_idx, int hap2_end_idx, 
           int start, int end,
           int phase, int orig_phase_dist, int swap_phase_dist);

    // pointers to variant/cluster data
    std::shared_ptr<variantCalls> cal1_vars = nullptr; 
    std::shared_ptr<variantCalls> cal2_vars = nullptr; 
    std::shared_ptr<variantCalls> hap1_vars = nullptr;
    std::shared_ptr<variantCalls> hap2_vars = nullptr;

    // cluster indices of superclusters
    std::vector<int> cal1_beg_idx, cal1_end_idx, cal2_beg_idx, cal2_end_idx,
        hap1_beg_idx, hap1_end_idx, hap2_beg_idx, hap2_end_idx;
    std::vector<int> starts, ends; // ref positions

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

void cluster(std::unique_ptr<vcfData> & vcf);

#endif
