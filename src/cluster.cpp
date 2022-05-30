#include <string>
#include <vector>

#include "cluster.h"
#include "vcf.h"

void ctgClusters::add(
           variantCalls* cal1_vars, int cal1_beg_idx, int cal1_end_idx, 
           variantCalls* cal2_vars, int cal2_beg_idx, int cal2_end_idx, 
           variantCalls* hap1_vars, int hap1_beg_idx, int hap1_end_idx, 
           variantCalls* hap2_vars, int hap2_beg_idx, int hap2_end_idx, 
           int phase, int orig_phase_dist, int swap_phase_dist) {

    // save pointers to variant data
    this->cal1_vars = cal1_vars;
    this->cal1_beg_idx.push_back(cal1_beg_idx);
    this->cal1_end_idx.push_back(cal1_end_idx);
    this->cal2_vars = cal1_vars;
    this->cal2_beg_idx.push_back(cal1_beg_idx);
    this->cal2_end_idx.push_back(cal1_end_idx);
    this->hap1_vars = hap1_vars;
    this->hap1_beg_idx.push_back(hap1_beg_idx);
    this->hap1_end_idx.push_back(hap1_end_idx);
    this->hap2_vars = hap1_vars;
    this->hap2_beg_idx.push_back(hap1_beg_idx);
    this->hap2_end_idx.push_back(hap1_end_idx);

    // store phasing info
    this->phase.push_back(phase);
    this->orig_phase_dist.push_back(orig_phase_dist);
    this->swap_phase_dist.push_back(swap_phase_dist);

    this->n++;
}
