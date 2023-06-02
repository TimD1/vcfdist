#include <string>
#include <vector>

#include "globals.h"
#include "dist.h"
#include "cluster.h"
#include "print.h"

ctgSuperclusters::ctgSuperclusters() {
    this->superclusters = std::vector< std::vector< std::vector<int> > > (2,
            std::vector< std::vector<int> > (2, std::vector<int>()));
    this->ctg_variants = std::vector< std::vector< std::shared_ptr<ctgVariants> > > (2, 
            std::vector< std::shared_ptr<ctgVariants> >(2, nullptr));
}

void ctgSuperclusters::add_supercluster(
           std::vector<int> brks, int beg, int end) {
    for (int i = 0; i < CALLSETS*HAPS; i++)
        this->superclusters[i>>1][i&1].push_back(brks[i]);
    this->begs.push_back(beg);
    this->ends.push_back(end);
    this->n++;
}

void ctgSuperclusters::add_phasing(
           int phase, 
           int orig_phase_dist, 
           int swap_phase_dist) {
    this->phase.push_back(phase);
    this->orig_phase_dist.push_back(orig_phase_dist);
    this->swap_phase_dist.push_back(swap_phase_dist);
}

/******************************************************************************/

superclusterData::superclusterData(
        std::unique_ptr<variantData> & query_ptr,
        std::unique_ptr<variantData> & truth_ptr,
        std::shared_ptr<fastaData> ref_ptr) {

    // set reference pointer
    this->ref = ref_ptr;

    // create list of all contigs covered by truth/query
    for (int i = 0; i < int(query_ptr->contigs.size()); i++) {
        std::string ctg = query_ptr->contigs[i];
        if (std::find(this->contigs.begin(), this->contigs.end(), ctg) == 
                this->contigs.end()) {
            this->contigs.push_back(ctg);
            this->lengths.push_back(query_ptr->lengths[i]);
            this->ploidy.push_back(query_ptr->ploidy[i]);
            this->ctg_superclusters[ctg] = std::shared_ptr<ctgSuperclusters>(
                    new ctgSuperclusters());
        }
    }
    for (int i = 0; i < int(truth_ptr->contigs.size()); i++) {
        std::string ctg = truth_ptr->contigs[i];
        if (std::find(this->contigs.begin(), this->contigs.end(), ctg) == 
                this->contigs.end()) {
            this->contigs.push_back(ctg);
            this->lengths.push_back(truth_ptr->lengths[i]);
            this->ploidy.push_back(truth_ptr->ploidy[i]);
            this->ctg_superclusters[ctg] = std::shared_ptr<ctgSuperclusters>(
                    new ctgSuperclusters());
        }
    }

    // set pointers to variant lists (per contig)
    for (std::string ctg : this->contigs) {
        try {
            this->ctg_superclusters[ctg]->ctg_variants[QUERY][HAP1] = 
                    query_ptr->ctg_variants[HAP1][ctg];
            this->ctg_superclusters[ctg]->ctg_variants[QUERY][HAP2] = 
                    query_ptr->ctg_variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Query VCF does not contain contig '%s'", ctg.data());
        }
        try {
            this->ctg_superclusters[ctg]->ctg_variants[TRUTH][HAP1] = 
                    truth_ptr->ctg_variants[HAP1][ctg];
            this->ctg_superclusters[ctg]->ctg_variants[TRUTH][HAP2] = 
                    truth_ptr->ctg_variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Truth VCF does not contain contig '%s'", ctg.data());
        }
    }

    // generate variant clusters
    this->gap_supercluster();
}


/******************************************************************************/


void superclusterData::gap_supercluster() {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("Superclustering truth and call variants across haplotypes");

    // iterate over each contig
    int total_superclusters = 0;
    int largest_supercluster = 0;
    int total_vars = 0;
    int most_vars = 0;
    int total_bases = 0;
    for (std::string ctg : this->contigs) {

        // skip empty contigs
        int nvars = 0;
        for (int i = 0; i < CALLSETS*HAPS; i++)
            nvars += this->ctg_superclusters[ctg]->ctg_variants[i>>1][i&1]->n;
        if (nvars > 0) { 
            if (g.verbosity >= 1) INFO("  Contig '%s': %d variants", ctg.data(), nvars); 
        } else { continue; }

        // for each cluster of variants (merge query and truth haps)
        auto vars = this->ctg_superclusters[ctg]->ctg_variants;
        std::vector<int> brks = {0, 0, 0, 0}; // start of current supercluster
        while (true) {

            // init: empty supercluster
            std::vector<int> next_brks = brks; // end of current supercluster
            std::vector<int> poss(4, std::numeric_limits<int>::max());
            for (int i = 0; i < CALLSETS*HAPS; i++) {
                if (brks[i] < int(vars[i>>1][i&1]->clusters.size())-1) {
                    poss[i] = vars[i>>1][i&1]->poss[ 
                        vars[i>>1][i&1]->clusters[ next_brks[i] ] ];
                }
            }

            // get first cluster, end if all haps are off end
            int idx = std::distance(poss.begin(),
                    std::min_element(poss.begin(), poss.end()));
            if (poss[idx] == std::numeric_limits<int>::max()) break;
            next_brks[idx]++;

            // initialize cluster merging with first to start
            int curr_end_pos = vars[idx>>1][idx&1]->poss[
                vars[idx>>1][idx&1]->clusters[next_brks[idx]]-1] +
                vars[idx>>1][idx&1]->rlens[
                    vars[idx>>1][idx&1]->clusters[next_brks[idx]]-1] + 1;
            poss[idx] = next_brks[idx] < int(vars[idx>>1][idx&1]->clusters.size())-1 ?
                vars[idx>>1][idx&1]->poss[
                    vars[idx>>1][idx&1]->clusters[next_brks[idx]]]-1 :
                    std::numeric_limits<int>::max();

            // keep expanding cluster while possible
            bool just_merged = true;
            while (just_merged) {
                just_merged = false;
                for (int i = 0; i < CALLSETS*HAPS; i++) {
                    while (poss[i] < curr_end_pos + g.cluster_min_gap) {
                        next_brks[i]++;
                        curr_end_pos = std::max(curr_end_pos,
                            vars[i>>1][i&1]->poss[
                                vars[i>>1][i&1]->clusters[next_brks[i]]-1] +
                            vars[i>>1][i&1]->rlens[
                                vars[i>>1][i&1]->clusters[next_brks[i]]-1] + 1);
                        poss[i] = next_brks[i] < int(vars[i>>1][i&1]->clusters.size())-1 ?
                            vars[i>>1][i&1]->poss[
                                vars[i>>1][i&1]->clusters[next_brks[i]]]-1 :
                                std::numeric_limits<int>::max();
                        just_merged = true;
                    }
                }
            }

            // get supercluster start/end positions (allowing empty haps)
            int beg_pos = std::numeric_limits<int>::max();
            int end_pos = -1;
            for (int i = 0; i < CALLSETS*HAPS; i++) {
                if (next_brks[i] - brks[i]) {
                    beg_pos = std::min(beg_pos,
                        vars[i>>1][i&1]->poss[
                            vars[i>>1][i&1]->clusters[brks[i]]]-1);
                    end_pos = std::max(end_pos,
                            vars[i>>1][i&1]->poss[
                                vars[i>>1][i&1]->clusters[next_brks[i]]-1] +
                            vars[i>>1][i&1]->rlens[
                                vars[i>>1][i&1]->clusters[next_brks[i]]-1] + 1);
                }
            }

            // update summary metrics
            largest_supercluster = std::max(largest_supercluster, end_pos-beg_pos);
            total_bases += end_pos-beg_pos;
            int this_vars = 0;
            for (int i = 0; i < CALLSETS*HAPS; i++) {
                this_vars += next_brks[i] - brks[i];
            }
            most_vars = std::max(most_vars, this_vars);
            total_vars += this_vars;

            // debug print
            if (false) {
                printf("\nSUPERCLUSTER: %d\n", this->ctg_superclusters[ctg]->n);
                printf("POS: %d-%d\n", beg_pos, end_pos);
                printf("SIZE: %d\n", end_pos - beg_pos);
                for (int i = 0; i < CALLSETS*HAPS; i++) {
                    printf("%s%d: clusters %d-%d, %d total\n",
                            callset_strs[i>>1].data(), (i&1)+1,
                            brks[i], next_brks[i], 
                            std::max(int(vars[i>>1][i&1]->clusters.size())-1, 0));
                }
            }

            // save alignment information
            this->ctg_superclusters[ctg]->add_supercluster(brks, beg_pos, end_pos);

            // reset for next merged cluster
            brks = next_brks;
        }

        // add sentinel, not actually a supercluster
        this->ctg_superclusters[ctg]->add_supercluster(brks, 
            std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
        this->ctg_superclusters[ctg]->n--;

        if (this->ctg_superclusters[ctg]->n && g.verbosity >= 1)
            INFO("    %d superclusters", this->ctg_superclusters[ctg]->n);
        total_superclusters += this->ctg_superclusters[ctg]->n;
    }
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("           Total superclusters: %d", total_superclusters);
    if (g.verbosity >= 1) INFO("  Largest supercluster (bases): %d", largest_supercluster);
    if (g.verbosity >= 1) INFO("  Largest supercluster  (vars): %d", most_vars);
    if (g.verbosity >= 1) INFO("  Average supercluster (bases): %d", total_bases / total_superclusters);
    if (g.verbosity >= 1) INFO("  Average supercluster  (vars): %d", total_vars / total_superclusters);
}


/******************************************************************************/

/* Cluster variants based on minimum gap length for independence */
void gap_cluster(std::unique_ptr<variantData> & vcf, int callset) {
    /* Add single-VCF cluster indices to `variantData` */
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("Gap Clustering %s VCF '%s'", 
            callset_strs[callset].data(), vcf->filename.data());

    // cluster each contig
    int largest_cluster_vars = 0;
    int largest_cluster_bases = 0;
    for (std::string ctg : vcf->contigs) {

        // only print for non-empty contigs
        if (vcf->ctg_variants[0][ctg]->n + vcf->ctg_variants[1][ctg]->n)
            if (g.verbosity >= 1) INFO("  Contig '%s'", ctg.data());

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < HAPS; hap++) {
            int nvar = vcf->ctg_variants[hap][ctg]->n;
            if (nvar && g.verbosity >= 1) INFO("    Haplotype %d", hap+1);
            int prev_end = -g.cluster_min_gap * 2;
            int cluster_start = 0;
            int pos = 0;
            int end = 0;
            int var_idx = 0;
            int prev_var_idx = 0;
            for (; var_idx < vcf->ctg_variants[hap][ctg]->n; var_idx++) {
                pos = vcf->ctg_variants[hap][ctg]->poss[var_idx];
                end = pos + vcf->ctg_variants[hap][ctg]->rlens[var_idx];
                if (pos - prev_end > g.cluster_min_gap) {
                    // update summary metrics
                    largest_cluster_bases = std::max(largest_cluster_bases, prev_end-cluster_start);
                    largest_cluster_vars = std::max(largest_cluster_vars, var_idx - prev_var_idx);
                    prev_var_idx = var_idx;
                    cluster_start = pos;
                    // add cluster
                    vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
                }
                prev_end = std::max(prev_end, end);
            }
            vcf->ctg_variants[hap][ctg]->add_cluster(var_idx);
            if (nvar && g.verbosity >= 1) INFO("      %d clusters", 
                    int(vcf->ctg_variants[hap][ctg]->clusters.size()-1));
        }
    }
    if (g.verbosity >= 1) INFO("  Largest cluster (bases): %d", largest_cluster_bases);
    if (g.verbosity >= 1) INFO("  Largest cluster  (vars): %d", largest_cluster_vars);
}


/******************************************************************************/


/* Add single-VCF cluster indices to `variantData`. This version assumes that
 * all variant calls are true positives (doesn't allow skipping)
 */
void swg_cluster(std::unique_ptr<variantData> & vcf, 
        int sub, int open, int extend, int callset, bool print /* = false */) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("Smith-Waterman Clustering %s VCF '%s'", 
            callset_strs[callset].data(), vcf->filename.data());

    // cluster each contig
    for (std::string ctg : vcf->contigs) {

        // only print for non-empty contigs
        if (vcf->ctg_variants[0][ctg]->n + vcf->ctg_variants[1][ctg]->n) {
            if (g.verbosity >= 1) INFO("  Contig '%s'", ctg.data());
        } else { continue; }

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < HAPS; hap++) {

            // init: each variant is its own cluster
            int nvar = vcf->ctg_variants[hap][ctg]->n;
            if (nvar) { 
                if (g.verbosity >= 1) INFO("    Haplotype %d", hap+1); 
            }
            else { continue; }
            std::vector<int> prev_clusters(nvar+1);
            for (int i = 0; i < nvar+1; i++) prev_clusters[i] = i;
            std::vector<bool> prev_merged(nvar+1, true);

            std::vector<int> right_reach, left_reach;
            std::vector<int> next_clusters;
            std::vector<bool> next_merged;

            // while clusters are being merged, loop
            int iter = 0;
            while (std::find(prev_merged.begin(), 
                        prev_merged.end(), true) != prev_merged.end()) {
                iter++;
                if (iter > g.max_cluster_itrs) break;
                if (print) printf("Iteration %d\n", iter);

                // count clusters currently being expanded
                int active = 0;
                for (size_t i = 0; i < prev_merged.size()-1; i++) {
                    if (prev_merged[i]) active++;
                }
                if (g.verbosity >= 2)
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

                    // no left cluster, not actual cluster
                    if (clust == 0 || clust == prev_clusters.size()-1) {
                        left_compute = false;
                    }
                    // no right cluster (second-to-last is actual last cluster)
                    if (clust >= prev_clusters.size()-2) {
                        right_compute = false;
                    }

                    // debug print
                    if (print && clust < prev_clusters.size()-1) {
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
                        score = calc_vcf_swg_score(
                                vcf->ctg_variants[hap][ctg], clust, clust+1,
                                sub, open, extend);
                        if (print) printf("orig score: %d\n", score);
                    }

                    // LEFT REACH 
                    // upper bound on distance off main diagonal
                    int buffer = score/extend + 1;

                    if (left_compute) { // calculate left reach
                        std::string query, ref;
                        std::vector< std::vector<int> > query_ref_ptrs, ref_query_ptrs;

                        // error checking
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust]-1 < 0)
                            ERROR("left var_idx < 0");
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust]-1 >= nvar)
                            ERROR("left var_idx >= nvar");
                        if (clust+1 >= vcf->ctg_variants[hap][ctg]->clusters.size())
                            ERROR("left next clust_idx >= nclust");
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1 < 0)
                            ERROR("left next var_idx < 0");
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1 >= nvar)
                            ERROR("left next var_idx >= nvar");

                        // just after last variant in this cluster
                        int end_pos = 
                            vcf->ctg_variants[hap][ctg]->poss[
                                vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1] +
                            vcf->ctg_variants[hap][ctg]->rlens[
                                vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1]+1;
                        // just after last variant in previous cluster
                        int beg_pos = std::max(0, 
                            vcf->ctg_variants[hap][ctg]->poss[
                                vcf->ctg_variants[hap][ctg]->clusters[clust]-1] +
                            vcf->ctg_variants[hap][ctg]->rlens[
                                vcf->ctg_variants[hap][ctg]->clusters[clust]-1] + 1
                            - buffer);

                        // generate reversed pointers/strings
                        generate_ptrs_strs(query, ref,
                                query_ref_ptrs, ref_query_ptrs, 
                                vcf->ctg_variants[hap][ctg],
                                clust, clust+1, beg_pos, end_pos,
                                vcf->ref, ctg 
                        );
                        reverse_ptrs_strs(query, ref, 
                                query_ref_ptrs, ref_query_ptrs);
                        
                        // get reference end pos of last variant in this cluster
                        int ref_section_start = end_pos - vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust]];

                        // calculate max reaching path to left
                        int reach = swg_max_reach(query, ref, 
                                query_ref_ptrs, ref_query_ptrs, 
                                sub, open, extend, score, 
                                print, true, ref_section_start); // reverse
                        l_reach = end_pos - reach;
                        if (print) printf("left reach: %d\n", reach);

                        if (print) {
                            printf("REF:        %s\n", ref.data());
                            printf("QUERY:      %s\n", query.data());
                            printf("ref start:  %d\n", ref_section_start);
                            printf("QUERY->REF: ");
                            for(size_t i = 0; i < query_ref_ptrs.size(); i++) 
                                printf("%d ", query_ref_ptrs[PTRS][i]); 
                            printf("\n");
                            printf("REF->QUERY: ");
                            for(size_t i = 0; i < ref_query_ptrs.size(); i++) 
                                printf("%d ", ref_query_ptrs[PTRS][i]); 
                            printf("\n");
                        }
                        
                    } else {
                        // past farthest right (unused)
                        l_reach = vcf->ctg_variants[hap][ctg]->poss[nvar-1]+10;
                    }
                    left_reach.push_back(l_reach);

                    // RIGHT REACH
                    if (right_compute) { // calculate right reach
                        std::string query, ref;
                        std::vector< std::vector<int> > query_ref_ptrs, ref_query_ptrs;

                        // error checking
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust] < 0)
                            ERROR("right var_idx < 0");
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust] >= nvar)
                            ERROR("right var_idx >= nvar");
                        if (clust+1 >= vcf->ctg_variants[hap][ctg]->clusters.size())
                            ERROR("right next clust_idx >= nclust");
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust+1] < 0)
                            ERROR("right next var_idx < 0");
                        if (vcf->ctg_variants[hap][ctg]->clusters[clust+1] >= nvar)
                            ERROR("right next var_idx >= nvar");

                        // right before current cluster
                        int beg_pos = vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust] ]-1;
                        // right before next cluster
                        int end_pos = vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]]
                                    + buffer;

                        generate_ptrs_strs(query, ref,
                                query_ref_ptrs, ref_query_ptrs, 
                                vcf->ctg_variants[hap][ctg], 
                                clust, clust+1, beg_pos, end_pos,
                                vcf->ref, ctg 
                        );

                        // get reference end pos of last variant in this cluster
                        int ref_section_start = vcf->ctg_variants[hap][ctg]->poss[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1]
                                    + vcf->ctg_variants[hap][ctg]->rlens[ 
                                    vcf->ctg_variants[hap][ctg]->clusters[clust+1]-1] - beg_pos;

                        // calculate max reaching path to right
                        int reach = swg_max_reach(query, ref, 
                                query_ref_ptrs, ref_query_ptrs, 
                                sub, open, extend, score, print, 
                                false, ref_section_start);
                        r_reach = beg_pos + reach;
                        if (print) printf("right reach: %d\n", reach);

                        if (print) {
                            printf("REF:        %s\n", ref.data());
                            printf("QUERY:      %s\n", query.data());
                            printf("ref_start:  %d\n", ref_section_start);
                            printf("QUERY->REF: ");
                            for(size_t i = 0; i < query_ref_ptrs.size(); i++) 
                                printf("%d ", query_ref_ptrs[PTRS][i]); 
                            printf("\n");
                            printf("REF->QUERY: ");
                            for(size_t i = 0; i < ref_query_ptrs.size(); i++) 
                                printf("%d ", ref_query_ptrs[PTRS][i]); 
                            printf("\n");
                        }
                        
                    } else { // non-adjacent, don't really compute
                        r_reach = -10; // past farthest left (unused)
                    }
                    right_reach.push_back(r_reach);
                    if (print) printf("span: %s - %s\n", 
                                l_reach == vcf->ctg_variants[hap][ctg]->poss[nvar-1]+10 ? 
                                    "X" : std::to_string(l_reach).data(), 
                                r_reach == -10 ? "X" : std::to_string(r_reach).data());
                }

                // merge dependent clusters
                size_t clust = 0;
                while (clust < prev_clusters.size()) {
                    int clust_size = 1;
                    int max_reach = right_reach[clust];
                    while (clust+clust_size < prev_clusters.size() &&  // merge
                            max_reach >= left_reach[clust+clust_size]) {
                        max_reach = std::max(max_reach,
                                right_reach[clust+clust_size]);
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

            if (nvar && g.verbosity >= 1) INFO("      %d resulting clusters.", int(prev_clusters.size()-1));
        }
    }
}
