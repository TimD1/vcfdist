#include <string>
#include <vector>

#include "globals.h"
#include "cluster.h"

void ctgClusters::set_variants(
           std::shared_ptr<ctgVariants> calls1_vars,
           std::shared_ptr<ctgVariants> calls2_vars,
           std::shared_ptr<ctgVariants> truth1_vars,
           std::shared_ptr<ctgVariants> truth2_vars) {
    this->calls1_vars = calls1_vars;
    this->calls2_vars = calls2_vars;
    this->truth1_vars = truth1_vars;
    this->truth2_vars = truth2_vars;
}

void ctgClusters::add_supercluster(
           int calls1_beg_idx, int calls1_end_idx, 
           int calls2_beg_idx, int calls2_end_idx, 
           int truth1_beg_idx, int truth1_end_idx, 
           int truth2_beg_idx, int truth2_end_idx, 
           int beg, int end,
           int phase, int orig_phase_dist, int swap_phase_dist) {

    // save pointers to variant data
    this->calls1_beg_idx.push_back(calls1_beg_idx);
    this->calls1_end_idx.push_back(calls1_end_idx);
    this->calls2_beg_idx.push_back(calls2_beg_idx);
    this->calls2_end_idx.push_back(calls2_end_idx);
    this->truth1_beg_idx.push_back(truth1_beg_idx);
    this->truth1_end_idx.push_back(truth1_end_idx);
    this->truth2_beg_idx.push_back(truth2_beg_idx);
    this->truth2_end_idx.push_back(truth2_end_idx);
    this->begs.push_back(beg);
    this->ends.push_back(end);

    // store phasing info
    this->phase.push_back(phase);
    this->orig_phase_dist.push_back(orig_phase_dist);
    this->swap_phase_dist.push_back(swap_phase_dist);

    this->n++;
}

/******************************************************************************/

void cluster(std::unique_ptr<variantData> & vcf) {
    /* Add single-VCF cluster indices to `variantData` */

    // cluster each contig
    for (std::string ctg : vcf->contigs) {

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < 2; hap++) {
            int prev_end = -g.gap * 2;
            int pos = 0;
            int end = 0;
            size_t var_idx = 0;
            for (; 
                    var_idx < vcf->ctg_variants[hap][ctg]->poss.size(); var_idx++) {
                pos = vcf->ctg_variants[hap][ctg]->poss[var_idx];
                end = pos + vcf->ctg_variants[hap][ctg]->rlens[var_idx];
                if (pos - prev_end > g.gap)
                    vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
                prev_end = std::max(prev_end, end);
            }
            vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
        }
    }
}
