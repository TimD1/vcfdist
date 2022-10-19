#include <string>
#include <vector>

#include "globals.h"
#include "cluster.h"
#include "print.h"
#include "dist.h"

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
           int beg, int end) {
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
    this->n++;
}

void ctgClusters::add_phasing(
           int phase, 
           int orig_phase_dist, 
           int swap_phase_dist) {
    this->phase.push_back(phase);
    this->orig_phase_dist.push_back(orig_phase_dist);
    this->swap_phase_dist.push_back(swap_phase_dist);
}

/******************************************************************************/

clusterData::clusterData(
        std::unique_ptr<variantData> & calls_ptr,
        std::unique_ptr<variantData> & truth_ptr,
        std::shared_ptr<fastaData> ref_ptr) {

    // set reference pointer
    this->ref = ref_ptr;

    // create list of all contigs covered by truth/calls
    for (std::string ctg : calls_ptr->contigs) {
        if (std::find(this->contigs.begin(), this->contigs.end(), ctg) == this->contigs.end()) {
            this->contigs.push_back(ctg);
            this->ctg_superclusters[ctg] = std::shared_ptr<ctgClusters>(new ctgClusters());
        }
    }
    for (std::string ctg : truth_ptr->contigs) {
        if (std::find(this->contigs.begin(), this->contigs.end(), ctg) == this->contigs.end()) {
            this->contigs.push_back(ctg);
            this->ctg_superclusters[ctg] = std::shared_ptr<ctgClusters>(new ctgClusters());
        }
    }

    // set pointers to variant lists (per contig)
    for (std::string ctg : this->contigs) {
        std::shared_ptr<ctgVariants> calls1_vars, calls2_vars, truth1_vars, truth2_vars;
        try {
            calls1_vars = calls_ptr->ctg_variants[HAP1][ctg];
            calls2_vars = calls_ptr->ctg_variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Calls VCF does not contain contig '%s'", ctg.data());
        }
        try {
            truth1_vars = truth_ptr->ctg_variants[HAP1][ctg];
            truth2_vars = truth_ptr->ctg_variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Truth VCF does not contain contig '%s'", ctg.data());
        }
        this->ctg_superclusters[ctg]->set_variants(
                calls1_vars, calls2_vars, truth1_vars, truth2_vars);
    }

    // generate variant clusters
    this->gap_supercluster();
}


/******************************************************************************/


void clusterData::gap_supercluster() {

    // iterate over each contig
    for (std::string ctg : this->contigs) {

        // set shorter var names
        std::shared_ptr<ctgVariants> calls1_vars = this->ctg_superclusters[ctg]->calls1_vars;
        std::shared_ptr<ctgVariants> calls2_vars = this->ctg_superclusters[ctg]->calls2_vars;
        std::shared_ptr<ctgVariants> truth1_vars = this->ctg_superclusters[ctg]->truth1_vars;
        std::shared_ptr<ctgVariants> truth2_vars = this->ctg_superclusters[ctg]->truth2_vars;

        // for each cluster of variants (merge calls and truth haps)
        size_t calls1_clust_beg_idx = 0;
        size_t calls2_clust_beg_idx = 0;
        size_t truth1_clust_beg_idx = 0;
        size_t truth2_clust_beg_idx = 0;
        while (calls1_clust_beg_idx < calls1_vars->clusters.size()-1 || // last cluster added for end
               calls2_clust_beg_idx < calls2_vars->clusters.size()-1 ||
               truth1_clust_beg_idx < truth1_vars->clusters.size()-1 ||
               truth2_clust_beg_idx < truth2_vars->clusters.size()-1) {

            // init: empty supercluster
            size_t calls1_clust_end_idx = calls1_clust_beg_idx;
            size_t calls2_clust_end_idx = calls2_clust_beg_idx;
            size_t truth1_clust_end_idx = truth1_clust_beg_idx;
            size_t truth2_clust_end_idx = truth2_clust_beg_idx;

            int calls1_pos = (calls1_clust_beg_idx < calls1_vars->clusters.size()-1) ? 
                    calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int calls2_pos = (calls2_clust_beg_idx < calls2_vars->clusters.size()-1) ? 
                    calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int truth1_pos = (truth1_clust_beg_idx < truth1_vars->clusters.size()-1) ? 
                    truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int truth2_pos = (truth2_clust_beg_idx < truth2_vars->clusters.size()-1) ? 
                    truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();

            // initialize cluster merging with first to start
            int curr_end_pos = 0;
            if (calls1_pos <= truth1_pos && calls1_pos <= calls2_pos && calls1_pos <= truth2_pos) {
                calls1_clust_end_idx += 1;
                curr_end_pos = calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]-1] +
                        calls1_vars->rlens[calls1_vars->clusters[calls1_clust_end_idx]-1] + 1;
                calls1_pos = (calls1_clust_end_idx < calls1_vars->clusters.size()-1) ? 
                    calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (calls2_pos <= truth1_pos && calls2_pos <= calls1_pos && 
                    calls2_pos <= truth2_pos) {
                calls2_clust_end_idx += 1;
                curr_end_pos = calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]-1] +
                        calls2_vars->rlens[calls2_vars->clusters[calls2_clust_end_idx]-1] + 1;
                calls2_pos = (calls2_clust_end_idx < calls2_vars->clusters.size()-1) ? 
                    calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (truth1_pos <= calls1_pos && truth1_pos <= calls2_pos && 
                    truth1_pos <= truth2_pos) {
                truth1_clust_end_idx += 1;
                curr_end_pos = truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]-1] +
                        truth1_vars->rlens[truth1_vars->clusters[truth1_clust_end_idx]-1] + 1;
                truth1_pos = (truth1_clust_end_idx < truth1_vars->clusters.size()-1) ? 
                    truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else {
                truth2_clust_end_idx += 1;
                curr_end_pos = truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] +
                        truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1] + 1;
                truth2_pos = (truth2_clust_end_idx < truth2_vars->clusters.size()-1) ? 
                    truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            }


            // keep expanding cluster while possible
            bool just_merged = true;
            while (just_merged) {
                just_merged = false;
                while (truth1_pos < curr_end_pos + g.cluster_min_gap) {
                    truth1_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]-1] + 
                            truth1_vars->rlens[truth1_vars->clusters[truth1_clust_end_idx]-1] + 1);
                    truth1_pos = (truth1_clust_end_idx < truth1_vars->clusters.size()-1) ? 
                        truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (truth2_pos < curr_end_pos + g.cluster_min_gap) {
                    truth2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] + 
                            truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1] + 1);
                    truth2_pos = (truth2_clust_end_idx < truth2_vars->clusters.size()-1) ? 
                        truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (calls1_pos < curr_end_pos + g.cluster_min_gap) {
                    calls1_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos, 
                            calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]-1] + 
                            calls1_vars->rlens[calls1_vars->clusters[calls1_clust_end_idx]-1] + 1);
                    calls1_pos = (calls1_clust_end_idx < calls1_vars->clusters.size()-1) ? 
                        calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (calls2_pos < curr_end_pos + g.cluster_min_gap) {
                    calls2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]-1] + 
                            calls2_vars->rlens[calls2_vars->clusters[calls2_clust_end_idx]-1] + 1);
                    calls2_pos = (calls2_clust_end_idx < calls2_vars->clusters.size()-1) ? 
                        calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
            }

            // get supercluster start/end positions (allowing empty haps)
            int beg_pos = std::numeric_limits<int>::max();
            int end_pos = -1;
            if (calls1_clust_end_idx - calls1_clust_beg_idx) { // calls1 vars present
                beg_pos = std::min(beg_pos, 
                        calls1_vars->poss[calls1_vars->clusters[calls1_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]-1] + 
                        calls1_vars->rlens[calls1_vars->clusters[calls1_clust_end_idx]-1]+1);
            }
            if (calls2_clust_end_idx - calls2_clust_beg_idx) { // calls2 vars present
                beg_pos = std::min(beg_pos, 
                        calls2_vars->poss[calls2_vars->clusters[calls2_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]-1] + 
                        calls2_vars->rlens[calls2_vars->clusters[calls2_clust_end_idx]-1]+1);
            }
            if (truth1_clust_end_idx - truth1_clust_beg_idx) { // truth1 vars present
                beg_pos = std::min(beg_pos, 
                        truth1_vars->poss[truth1_vars->clusters[truth1_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]-1] + 
                        truth1_vars->rlens[truth1_vars->clusters[truth1_clust_end_idx]-1]+1);
            }
            if (truth2_clust_end_idx - truth2_clust_beg_idx) { // truth2 vars present
                beg_pos = std::min(beg_pos, 
                        truth2_vars->poss[truth2_vars->clusters[truth2_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] + 
                        truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1]+1);
            }

            // save alignment information
            this->ctg_superclusters[ctg]->add_supercluster(
                    calls1_clust_beg_idx, calls1_clust_end_idx,
                    calls2_clust_beg_idx, calls2_clust_end_idx,
                    truth1_clust_beg_idx, truth1_clust_end_idx,
                    truth2_clust_beg_idx, truth2_clust_end_idx,
                    beg_pos, end_pos);

            // reset for next merged cluster
            calls1_clust_beg_idx = calls1_clust_end_idx;
            calls2_clust_beg_idx = calls2_clust_end_idx;
            truth1_clust_beg_idx = truth1_clust_end_idx;
            truth2_clust_beg_idx = truth2_clust_end_idx;
        }
    }
}


/******************************************************************************/


void cluster(std::unique_ptr<variantData> & vcf) {
    /* Add single-VCF cluster indices to `variantData` */

    // cluster each contig
    for (std::string ctg : vcf->contigs) {

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < 2; hap++) {
            int prev_end = -g.cluster_min_gap * 2;
            int pos = 0;
            int end = 0;
            size_t var_idx = 0;
            for (; 
                    var_idx < vcf->ctg_variants[hap][ctg]->poss.size(); var_idx++) {
                pos = vcf->ctg_variants[hap][ctg]->poss[var_idx];
                end = pos + vcf->ctg_variants[hap][ctg]->rlens[var_idx];
                if (pos - prev_end > g.cluster_min_gap)
                    vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
                prev_end = std::max(prev_end, end);
            }
            vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
        }
    }
}


/******************************************************************************/


/* Add single-VCF cluster indices to `variantData` */
void sw_cluster(std::unique_ptr<variantData> & vcf) {

    // cluster each contig
    for (std::string ctg : vcf->contigs) {

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < 2; hap++) {

            // init: each variant is its own cluster
            size_t nvar = vcf->ctg_variants[hap][ctg]->n;
            std::vector<int> prev_clusters(nvar+1);
            for (size_t i = 0; i < nvar+1; i++) prev_clusters[i] = i;
            std::vector<bool> prev_merged(nvar+1, true);
            prev_merged[nvar] = false; // end pointer isn't cluster

            std::vector<int> right_reach, left_reach;
            std::vector<int> next_clusters;
            std::vector<bool> next_merged;

            // while clusters are being merged, loop
            while (std::find(prev_merged.begin(), 
                        prev_merged.end(), true) != prev_merged.end()) {

                // save temp clustering (generate_ptrs_strs assumes clustered)
                vcf->ctg_variants[hap][ctg]->clusters = prev_clusters;

                // update all cluster reaches
                for (size_t clust = 0; clust < prev_clusters.size(); clust++) {

                    // only compute necessary reaches (adjacent merge)
                    bool left_compute = prev_merged[clust];
                    if (clust > 0)    
                        left_compute = left_compute || prev_merged[clust-1];
                    // can't merge first cluster left, last cluster is just end ptr
                    if (clust == 0 || clust >= prev_clusters.size()-1)
                        left_compute = false;
                    bool right_compute = prev_merged[clust];
                    if (clust < prev_merged.size()-1) 
                        right_compute = right_compute || prev_merged[clust+1];
                    // can't merge second to last cluster right, 
                    // last cluster is just end ptr
                    if (clust >= prev_clusters.size()-2)
                        right_compute = false;

                    // LEFT REACH
                    if (left_compute) {

                        // TODO: reverse align
                        
                    } else {
                        left_reach.push_back( // past farthest right
                                vcf->ctg_variants[hap][ctg]->poss[nvar-1]+1);
                    }

                    // RIGHT REACH
                    if (right_compute) {

                        printf("===========================================\n");
                        for (int var_idx = vcf->ctg_variants[hap][ctg]->clusters[clust]; 
                                var_idx < vcf->ctg_variants[hap][ctg]->clusters[clust+1]; var_idx++) {
                            printf("%s %d %s %s var:%d\n",
                                    ctg.data(), 
                                    vcf->ctg_variants[hap][ctg]->poss[var_idx],
                                    vcf->ctg_variants[hap][ctg]->refs[var_idx].data(),
                                    vcf->ctg_variants[hap][ctg]->alts[var_idx].data(),
                                    var_idx
                            );
                        }

                        // calculate right reach
                        std::string calls, ref, calls_str, ref_str;
                        std::vector<int> calls_ref_ptrs, ref_calls_ptrs;
                        generate_ptrs_strs(
                                calls, ref, calls_str, ref_str,
                                calls_ref_ptrs, ref_calls_ptrs, 
                                vcf->ctg_variants[hap][ctg], 
                                vcf->ctg_variants[hap][ctg],
                                clust, 0, clust+1, 0,
                                vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust]],
                                vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]],
                                vcf->ref, ctg 
                        );
                        
                        // calculate existing VCF score
                        int score = calc_vcf_sw_score(
                                vcf->ctg_variants[hap][ctg], clust, clust+1);
                        printf("orig score: %d\n", score);

                        /* // calculate max reaching path to right */
                        /* int reach = max_reach(calls, ref, */ 
                        /*         calls_ref_ptrs, ref_calls_ptrs, score); */
                        /* printf("reach: %d\n", reach); */

                        printf("clusters %d-%d, pos %d-%d, ctg %s nvar %d, nclust %d\n",
                                int(clust), int(clust+1), 
                                vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust]],
                                vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]],
                                ctg.data(), int(nvar), int(prev_clusters.size())
                                );

                        printf("REF:        %s\n", ref_str.data());
                        printf("CALLS:      %s\n", calls_str.data());
                        /* printf("CALLS->REF: "); */
                        /* for(size_t i = 0; i < calls_ref_ptrs.size(); i++) */ 
                        /*     printf("%d ", calls_ref_ptrs[i]); */ 
                        /* printf("\n"); */
                        /* printf("REF->CALLS: "); */
                        /* for(size_t i = 0; i < ref_calls_ptrs.size(); i++) */ 
                        /*     printf("%d ", ref_calls_ptrs[i]); */ 
                        /* printf("\n"); */
                        
                    } else { // non-adjacent, don't really compute
                        right_reach.push_back(-1); // past farthest left
                    }
                }
                ERROR("Breakpoint");

                // merge dependent clusters
                size_t clust = 0;
                while (clust < prev_clusters.size()) {
                    int next = 1;
                    while (clust+next < prev_clusters.size() && 
                            right_reach[clust] >= left_reach[clust+next]) { // merge
                        next++;
                    }
                    next_clusters.push_back(clust);
                    next_merged.push_back(next > 1); // true if merge occurred
                    clust += next;
                }

                // reset for next iteration
                prev_clusters = next_clusters;
                prev_merged = next_merged;
                next_clusters.clear();
                next_merged.clear();
                left_reach.clear();
                right_reach.clear();
            }

            // save final clustering
            vcf->ctg_variants[hap][ctg]->clusters = prev_clusters;
        }
    }
}
