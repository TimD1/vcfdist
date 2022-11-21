#include <string>
#include <vector>

#include "globals.h"
#include "cluster.h"
#include "print.h"
#include "dist.h"

void ctgClusters::add_supercluster(
           int query1_beg_idx, int query1_end_idx, 
           int query2_beg_idx, int query2_end_idx, 
           int truth1_beg_idx, int truth1_end_idx, 
           int truth2_beg_idx, int truth2_end_idx, 
           int beg, int end) {
    this->query1_beg_idx.push_back(query1_beg_idx);
    this->query1_end_idx.push_back(query1_end_idx);
    this->query2_beg_idx.push_back(query2_beg_idx);
    this->query2_end_idx.push_back(query2_end_idx);
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
        std::unique_ptr<variantData> & query_ptr,
        std::unique_ptr<variantData> & truth_ptr,
        std::shared_ptr<fastaData> ref_ptr) {

    // set reference pointer
    this->ref = ref_ptr;

    // create list of all contigs covered by truth/query
    for (std::string ctg : query_ptr->contigs) {
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
        try {
            this->ctg_superclusters[ctg]->query1_vars = query_ptr->ctg_variants[HAP1][ctg];
            this->ctg_superclusters[ctg]->query2_vars = query_ptr->ctg_variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Query VCF does not contain contig '%s'", ctg.data());
        }
        try {
            this->ctg_superclusters[ctg]->truth1_vars = truth_ptr->ctg_variants[HAP1][ctg];
            this->ctg_superclusters[ctg]->truth2_vars = truth_ptr->ctg_variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Truth VCF does not contain contig '%s'", ctg.data());
        }
    }

    // generate variant clusters
    this->gap_supercluster();
}


/******************************************************************************/


void clusterData::gap_supercluster() {
    INFO(" ");
    INFO("Superclustering truth and call variants across haplotypes");

    // iterate over each contig
    int total_superclusters = 0;
    int largest_supercluster = 0;
    for (std::string ctg : this->contigs) {

        // set shorter var names
        std::shared_ptr<ctgVariants> query1_vars = this->ctg_superclusters[ctg]->query1_vars;
        std::shared_ptr<ctgVariants> query2_vars = this->ctg_superclusters[ctg]->query2_vars;
        std::shared_ptr<ctgVariants> truth1_vars = this->ctg_superclusters[ctg]->truth1_vars;
        std::shared_ptr<ctgVariants> truth2_vars = this->ctg_superclusters[ctg]->truth2_vars;

        if (query1_vars->n + query2_vars->n + truth1_vars->n + truth2_vars->n) {
            INFO("  Contig '%s'", ctg.data());
        } else {continue;}

        // for each cluster of variants (merge query and truth haps)
        int query1_clust_beg_idx = 0;
        int query2_clust_beg_idx = 0;
        int truth1_clust_beg_idx = 0;
        int truth2_clust_beg_idx = 0;
        while (query1_clust_beg_idx < int(query1_vars->clusters.size())-1 || // last cluster added for end
               query2_clust_beg_idx < int(query2_vars->clusters.size())-1 ||
               truth1_clust_beg_idx < int(truth1_vars->clusters.size())-1 ||
               truth2_clust_beg_idx < int(truth2_vars->clusters.size())-1) {

            // init: empty supercluster
            int query1_clust_end_idx = query1_clust_beg_idx;
            int query2_clust_end_idx = query2_clust_beg_idx;
            int truth1_clust_end_idx = truth1_clust_beg_idx;
            int truth2_clust_end_idx = truth2_clust_beg_idx;

            int query1_pos = (query1_clust_beg_idx < int(query1_vars->clusters.size())-1) ? 
                    query1_vars->poss[query1_vars->clusters[query1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int query2_pos = (query2_clust_beg_idx < int(query2_vars->clusters.size())-1) ? 
                    query2_vars->poss[query2_vars->clusters[query2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int truth1_pos = (truth1_clust_beg_idx < int(truth1_vars->clusters.size())-1) ? 
                    truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int truth2_pos = (truth2_clust_beg_idx < int(truth2_vars->clusters.size())-1) ? 
                    truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();

            // initialize cluster merging with first to start
            int curr_end_pos = 0;
            if (query1_pos <= truth1_pos && query1_pos <= query2_pos && query1_pos <= truth2_pos) {
                query1_clust_end_idx += 1;
                curr_end_pos = query1_vars->poss[query1_vars->clusters[query1_clust_end_idx]-1] +
                        query1_vars->rlens[query1_vars->clusters[query1_clust_end_idx]-1] + 1;
                query1_pos = (query1_clust_end_idx < int(query1_vars->clusters.size())-1) ? 
                    query1_vars->poss[query1_vars->clusters[query1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (query2_pos <= truth1_pos && query2_pos <= query1_pos && 
                    query2_pos <= truth2_pos) {
                query2_clust_end_idx += 1;
                curr_end_pos = query2_vars->poss[query2_vars->clusters[query2_clust_end_idx]-1] +
                        query2_vars->rlens[query2_vars->clusters[query2_clust_end_idx]-1] + 1;
                query2_pos = (query2_clust_end_idx < int(query2_vars->clusters.size())-1) ? 
                    query2_vars->poss[query2_vars->clusters[query2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (truth1_pos <= query1_pos && truth1_pos <= query2_pos && 
                    truth1_pos <= truth2_pos) {
                truth1_clust_end_idx += 1;
                curr_end_pos = truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]-1] +
                        truth1_vars->rlens[truth1_vars->clusters[truth1_clust_end_idx]-1] + 1;
                truth1_pos = (truth1_clust_end_idx < int(truth1_vars->clusters.size())-1) ? 
                    truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else {
                truth2_clust_end_idx += 1;
                curr_end_pos = truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] +
                        truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1] + 1;
                truth2_pos = (truth2_clust_end_idx < int(truth2_vars->clusters.size())-1) ? 
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
                    truth1_pos = (truth1_clust_end_idx < int(truth1_vars->clusters.size())-1) ? 
                        truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (truth2_pos < curr_end_pos + g.cluster_min_gap) {
                    truth2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] + 
                            truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1] + 1);
                    truth2_pos = (truth2_clust_end_idx < int(truth2_vars->clusters.size())-1) ? 
                        truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (query1_pos < curr_end_pos + g.cluster_min_gap) {
                    query1_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos, 
                            query1_vars->poss[query1_vars->clusters[query1_clust_end_idx]-1] + 
                            query1_vars->rlens[query1_vars->clusters[query1_clust_end_idx]-1] + 1);
                    query1_pos = (query1_clust_end_idx < int(query1_vars->clusters.size())-1) ? 
                        query1_vars->poss[query1_vars->clusters[query1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (query2_pos < curr_end_pos + g.cluster_min_gap) {
                    query2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            query2_vars->poss[query2_vars->clusters[query2_clust_end_idx]-1] + 
                            query2_vars->rlens[query2_vars->clusters[query2_clust_end_idx]-1] + 1);
                    query2_pos = (query2_clust_end_idx < int(query2_vars->clusters.size())-1) ? 
                        query2_vars->poss[query2_vars->clusters[query2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
            }

            // get supercluster start/end positions (allowing empty haps)
            int beg_pos = std::numeric_limits<int>::max();
            int end_pos = -1;
            if (query1_clust_end_idx - query1_clust_beg_idx) { // query1 vars present
                beg_pos = std::min(beg_pos, 
                        query1_vars->poss[query1_vars->clusters[query1_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        query1_vars->poss[query1_vars->clusters[query1_clust_end_idx]-1] + 
                        query1_vars->rlens[query1_vars->clusters[query1_clust_end_idx]-1]+1);
            }
            if (query2_clust_end_idx - query2_clust_beg_idx) { // query2 vars present
                beg_pos = std::min(beg_pos, 
                        query2_vars->poss[query2_vars->clusters[query2_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        query2_vars->poss[query2_vars->clusters[query2_clust_end_idx]-1] + 
                        query2_vars->rlens[query2_vars->clusters[query2_clust_end_idx]-1]+1);
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
            largest_supercluster = std::max(largest_supercluster, end_pos-beg_pos);

            /* printf("\nPOS: %d-%d\n", beg_pos, end_pos); */
            /* printf("SIZE: %d\n", end_pos - beg_pos); */
            /* printf("QUERY1: clusters %d-%d of %d\n", query1_clust_beg_idx, */ 
            /*         query1_clust_end_idx, int(query1_vars->clusters.size())); */
            /* for(int i = query1_clust_beg_idx; i < query1_clust_end_idx; i++) */
            /*     printf("  %d", query1_vars->poss[query1_vars->clusters[i]); */
            /* printf("QUERY2: clusters %d-%d of %d\n", query2_clust_beg_idx, */ 
            /*         query2_clust_end_idx, int(query2_vars->clusters.size())); */
            /* printf("TRUTH1: clusters %d-%d of %d\n", truth1_clust_beg_idx, */ 
            /*         truth1_clust_end_idx, int(truth1_vars->clusters.size())); */
            /* printf("TRUTH2: clusters %d-%d of %d\n", truth2_clust_beg_idx, */ 
            /*         truth2_clust_end_idx, int(truth2_vars->clusters.size())); */

            // save alignment information
            this->ctg_superclusters[ctg]->add_supercluster(
                    query1_clust_beg_idx, query1_clust_end_idx,
                    query2_clust_beg_idx, query2_clust_end_idx,
                    truth1_clust_beg_idx, truth1_clust_end_idx,
                    truth2_clust_beg_idx, truth2_clust_end_idx,
                    beg_pos, end_pos);

            // reset for next merged cluster
            query1_clust_beg_idx = query1_clust_end_idx;
            query2_clust_beg_idx = query2_clust_end_idx;
            truth1_clust_beg_idx = truth1_clust_end_idx;
            truth2_clust_beg_idx = truth2_clust_end_idx;
        }

        if (this->ctg_superclusters[ctg]->n)
            INFO("    %d superclusters", this->ctg_superclusters[ctg]->n);
        total_superclusters += this->ctg_superclusters[ctg]->n;
    }
    INFO(" ");
    INFO("  Total superclusters:  %d", total_superclusters);
    INFO("  Largest supercluster: %d", largest_supercluster);
}


/******************************************************************************/


void cluster(std::unique_ptr<variantData> & vcf) {
    /* Add single-VCF cluster indices to `variantData` */
    INFO(" ");
    INFO("Gap Clustering VCF '%s'", vcf->filename.data());

    // cluster each contig
    for (std::string ctg : vcf->contigs) {

        // only print for non-empty contigs
        if (vcf->ctg_variants[0][ctg]->n + vcf->ctg_variants[1][ctg]->n)
            INFO("  Contig '%s'", ctg.data());

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < 2; hap++) {
            int nvar = vcf->ctg_variants[hap][ctg]->n;
            if (nvar) INFO("    Haplotype %d", hap+1);
            int prev_end = -g.cluster_min_gap * 2;
            int pos = 0;
            int end = 0;
            int var_idx = 0;
            for (; var_idx < vcf->ctg_variants[hap][ctg]->n; var_idx++) {
                pos = vcf->ctg_variants[hap][ctg]->poss[var_idx];
                end = pos + vcf->ctg_variants[hap][ctg]->rlens[var_idx];
                if (pos - prev_end > g.cluster_min_gap)
                    vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
                prev_end = std::max(prev_end, end);
            }
            vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
            if (nvar) INFO("      %d resulting clusters.", 
                    int(vcf->ctg_variants[hap][ctg]->clusters.size()-1));
        }
    }
}


/******************************************************************************/


/* Add single-VCF cluster indices to `variantData`. This version assumes that
 * all variant calls are true positives (doesn't allow skipping)
 */
void sw_cluster(std::unique_ptr<variantData> & vcf, int sub, int open, int extend) {
    INFO(" ");
    INFO("SW Clustering VCF '%s'", vcf->filename.data());

    // cluster each contig
    for (std::string ctg : vcf->contigs) {

        // only print for non-empty contigs
        if (vcf->ctg_variants[0][ctg]->n + vcf->ctg_variants[1][ctg]->n) {
            INFO("  Contig '%s'", ctg.data());
        } else { continue; }

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < 2; hap++) {

            // init: each variant is its own cluster
            size_t nvar = vcf->ctg_variants[hap][ctg]->n;
            if (nvar) INFO("    Haplotype %d", hap+1);
            std::vector<int> prev_clusters(nvar+1);
            for (size_t i = 0; i < nvar+1; i++) prev_clusters[i] = i;
            std::vector<bool> prev_merged(nvar+1, true);

            std::vector<int> right_reach, left_reach;
            std::vector<int> next_clusters;
            std::vector<bool> next_merged;

            // while clusters are being merged, loop
            int iter = 0;
            while (std::find(prev_merged.begin(), 
                        prev_merged.end(), true) != prev_merged.end()) {
                iter++;
                if (g.print_verbosity >= 2)
                    printf("Iteration %d\n", iter);

                // count clusters currently being expanded
                int active = 0;
                for (size_t i = 0; i < prev_merged.size()-1; i++) {
                    if (prev_merged[i]) active++;
                }
                INFO("      Iteration %d: %d clusters, %d active",
                        iter, int(prev_clusters.size()-1), active);

                // save temp clustering (generate_ptrs_strs assumes clustered)
                vcf->ctg_variants[hap][ctg]->clusters = prev_clusters;

                // update all cluster reaches
                for (size_t clust = 0; clust < prev_clusters.size(); clust++) {

                    // only compute necessary reaches (adjacent merge)
                    bool left_compute = prev_merged[clust];
                    if (clust > 0) // protect OOB
                        left_compute = left_compute || prev_merged[clust-1];
                    bool right_compute = prev_merged[clust];
                    if (clust < prev_merged.size() - 1) // protect OOB
                        right_compute = right_compute || prev_merged[clust+1];

                    // edge cases
                    if (clust == 0 || clust == prev_clusters.size()-1) {
                        left_compute = false;
                        right_compute = false;
                    } else if (clust == prev_clusters.size()-2) {
                        right_compute = false;
                    }

                    // debug print
                    if (g.print_verbosity >= 2 && clust < prev_clusters.size()-1) {
                        printf("\ncluster %d: vars %d-%d, pos %d-%d\n",
                                int(clust), 
                                vcf->ctg_variants[hap][ctg]->clusters[clust],
                                vcf->ctg_variants[hap][ctg]->clusters[clust+1],
                                vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust]],
                                vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1] +
                                vcf->ctg_variants[hap][ctg]->rlens[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1]);
                        for (int var_idx = vcf->ctg_variants[hap][ctg]->clusters[clust]; 
                                var_idx < vcf->ctg_variants[hap][ctg]->clusters[clust+1]; var_idx++) {
                            printf("    %s %d %s %s var:%d\n",
                                    ctg.data(), 
                                    vcf->ctg_variants[hap][ctg]->poss[var_idx],
                                    vcf->ctg_variants[hap][ctg]->refs[var_idx].size() ? 
                                        vcf->ctg_variants[hap][ctg]->refs[var_idx].data() : "_",
                                    vcf->ctg_variants[hap][ctg]->alts[var_idx].size() ? 
                                        vcf->ctg_variants[hap][ctg]->alts[var_idx].data() : "_",
                                    var_idx
                            );
                        }
                    }

                    int l_reach = 0, r_reach = 0, score = 0;
                    if (left_compute || right_compute) {
                        score = calc_vcf_sw_score(
                                vcf->ctg_variants[hap][ctg], clust, clust+1,
                                sub, open, extend);
                        /* printf("orig score: %d\n", score); */
                    }

                    // LEFT REACH

                    // TODO: remove buffer
                    // NOTE: this buffer is only temporary, to allow extension 
                    // past the next variant for cluster merging. It should 
                    // eventually be replaced with Djikstra->WFA max row.
                    int buffer = 10;
                    if (left_compute) {

                        // calculate left reach
                        std::string query, ref;
                        std::vector<int> query_ref_ptrs, ref_query_ptrs;
                        // just after last variant in previous cluster
                        int beg_pos = vcf->ctg_variants[hap][ctg]->poss[
                                vcf->ctg_variants[hap][ctg]->clusters[clust]-1] +
                            vcf->ctg_variants[hap][ctg]->rlens[
                                vcf->ctg_variants[hap][ctg]->clusters[clust]-1]+1
                                - buffer; // TODO: remove
                        // just after last variant in this cluster
                        int end_pos = vcf->ctg_variants[hap][ctg]->poss[
                                vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1] +
                            vcf->ctg_variants[hap][ctg]->rlens[
                                vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1]+1;
                        generate_ptrs_strs(query, ref,
                                query_ref_ptrs, ref_query_ptrs, 
                                vcf->ctg_variants[hap][ctg], 
                                vcf->ctg_variants[hap][ctg],
                                clust, 0, clust+1, 0, beg_pos, end_pos,
                                vcf->ref, ctg 
                        );
                        reverse_ptrs_strs(query, ref, 
                                query_ref_ptrs, ref_query_ptrs);
                        
                        // calculate max reaching path to left
                        int reach = sw_max_reach(query, ref, 
                                query_ref_ptrs, ref_query_ptrs, 
                                sub, open, extend,
                                score, true); // reverse
                        l_reach = end_pos - reach;
                        /* printf("left reach: %d\n", reach); */

                        /* printf("REF:        %s\n", ref.data()); */
                        /* printf("QUERY:      %s\n", query.data()); */
                        /* printf("QUERY->REF: "); */
                        /* for(size_t i = 0; i < query_ref_ptrs.size(); i++) */ 
                        /*     printf("%d ", query_ref_ptrs[i]); */ 
                        /* printf("\n"); */
                        /* printf("REF->QUERY: "); */
                        /* for(size_t i = 0; i < ref_query_ptrs.size(); i++) */ 
                        /*     printf("%d ", ref_query_ptrs[i]); */ 
                        /* printf("\n"); */
                        
                    } else {
                        // past farthest right (unused)
                        l_reach = vcf->ctg_variants[hap][ctg]->poss[nvar-1]+10;
                    }
                    left_reach.push_back(l_reach);

                    // RIGHT REACH
                    if (right_compute) {

                        // calculate right reach
                        std::string query, ref;
                        std::vector<int> query_ref_ptrs, ref_query_ptrs;
                        // right before current cluster
                        int beg_pos = vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust]]-1;
                        // right before next cluster
                        int end_pos = vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]]
                                    + buffer; // TODO: remove
                        generate_ptrs_strs(query, ref,
                                query_ref_ptrs, ref_query_ptrs, 
                                vcf->ctg_variants[hap][ctg], 
                                vcf->ctg_variants[hap][ctg],
                                clust, 0, clust+1, 0, beg_pos, end_pos,
                                vcf->ref, ctg 
                        );

                        // calculate max reaching path to right
                        int reach = sw_max_reach(query, ref, 
                                query_ref_ptrs, ref_query_ptrs, 
                                sub, open, extend, score);
                        r_reach = beg_pos + reach;
                        /* printf("right reach: %d\n", reach); */

                        /* printf("REF:        %s\n", ref.data()); */
                        /* printf("QUERY:      %s\n", query.data()); */
                        /* printf("QUERY->REF: "); */
                        /* for(size_t i = 0; i < query_ref_ptrs.size(); i++) */ 
                        /*     printf("%d ", query_ref_ptrs[i]); */ 
                        /* printf("\n"); */
                        /* printf("REF->QUERY: "); */
                        /* for(size_t i = 0; i < ref_query_ptrs.size(); i++) */ 
                        /*     printf("%d ", ref_query_ptrs[i]); */ 
                        /* printf("\n"); */
                        
                    } else { // non-adjacent, don't really compute
                        r_reach = -10; // past farthest left (unused)
                    }
                    right_reach.push_back(r_reach);
                    /* printf("span: %d-%d\n", l_reach, r_reach); */
                }

                // merge dependent clusters
                size_t clust = 0;
                while (clust < prev_clusters.size()) {
                    int clust_size = 1;
                    while (clust+clust_size < prev_clusters.size() &&  // merge
                            right_reach[clust] >= left_reach[clust+clust_size]) {
                        clust_size++;
                    }
                    next_clusters.push_back(prev_clusters[clust]);
                    next_merged.push_back(clust_size > 1); // true if merge occurred
                    clust += clust_size;
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

            if (nvar) INFO("      %d resulting clusters.", int(prev_clusters.size()-1));
        }
    }
}
