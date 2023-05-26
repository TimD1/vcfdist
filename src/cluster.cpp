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
            bool just_active = true;
            while (just_active) {
                just_active = false;
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
                        just_active = true;
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

            // reset for next active cluster
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
    if (g.verbosity >= 1 && total_superclusters) INFO("  Average supercluster (bases): %d", total_bases / total_superclusters);
    if (g.verbosity >= 1 && total_superclusters) INFO("  Average supercluster  (vars): %d", total_vars / total_superclusters);
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
void wf_swg_cluster(std::unique_ptr<variantData> & vcf, 
        int sub, int open, int extend, int callset, bool print /* = false */) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("Smith-Waterman Clustering %s VCF '%s'", 
            callset_strs[callset].data(), vcf->filename.data());

    // allocate this memory once, use on each cluster
    std::vector<int> offs_buffer(MATS * g.max_size*2 * std::max(sub, open+extend), -2);

    // cluster each contig
    for (std::string ctg : vcf->contigs) {

        // only print for non-empty contigs
        if (vcf->ctg_variants[0][ctg]->n + vcf->ctg_variants[1][ctg]->n) {
            if (g.verbosity >= 1) INFO("  Contig '%s'", ctg.data());
        } else { continue; }
        int len = vcf->lengths[ std::find(vcf->contigs.begin(), 
                vcf->contigs.end(), ctg) - vcf->contigs.begin() ];

        // cluster per-haplotype variants: vcf->ctg_variants[hap]
        for (int hap = 0; hap < HAPS; hap++) {

            // init: each variant is its own cluster
            auto vars = vcf->ctg_variants[hap][ctg];
            int nvar = vars->n;
            if (nvar) { 
                if (g.verbosity >= 1) INFO("    Haplotype %d", hap+1); 
            }
            else { continue; }

            // mark variant boundary between clusters (hence n+1), 
            // and whether it was active on the previous iteration.
            // initially, one variant per cluster
            std::vector<int> prev_clusters(nvar+1);
            for (int i = 0; i < nvar+1; i++) prev_clusters[i] = i;
            std::vector<bool> prev_active(nvar+1, true);

            std::vector<int> right_reach, left_reach;
            std::vector<int> next_clusters, tmp_clusters;
            std::vector<bool> next_active, tmp_active;

            // while clusters are being merged, loop
            int iter = 0;
            while (std::find(prev_active.begin(), 
                        prev_active.end(), true) != prev_active.end()) {
                iter++;
                if (iter > g.max_cluster_itrs) break;
                if (print) printf("Iteration %d\n", iter);

                // count clusters currently being expanded
                int active = 0;
                for (size_t i = 0; i < prev_active.size()-1; i++) {
                    if (prev_active[i]) active++;
                }
                if (g.verbosity >= 2)
                    INFO("      Iteration %d: %d clusters, %d active",
                        iter, int(prev_clusters.size()-1), active);

                // save temp clustering (generate_str assumes clustered)
                vars->clusters = prev_clusters;

                // update all cluster reaches
                for (size_t clust = 0; clust < prev_clusters.size(); clust++) {

                    // only compute necessary reaches (adjacent merge)
                    bool left_compute = prev_active[clust];
                    if (clust > 0) // protect OOB
                        left_compute = left_compute || prev_active[clust-1];
                    bool right_compute = prev_active[clust];
                    if (clust < prev_active.size() - 1) // protect OOB
                        right_compute = right_compute || prev_active[clust+1];

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
                        printf("\n\ncluster %d: vars %d-%d, pos %d-%d\n",
                                int(clust), vars->clusters[clust],
                                vars->clusters[clust+1],
                                vars->poss[vars->clusters[clust]],
                                vars->poss[vars->clusters[clust+1]-1] +
                                vars->rlens[vars->clusters[clust+1]-1]);
                        for (int var_idx = vars->clusters[clust]; 
                                var_idx < vars->clusters[clust+1]; var_idx++) {
                            printf("    %s %d %s %s var:%d\n",
                                    ctg.data(), 
                                    vars->poss[var_idx],
                                    vars->refs[var_idx].size() ? 
                                        vars->refs[var_idx].data() : "_",
                                    vars->alts[var_idx].size() ? 
                                        vars->alts[var_idx].data() : "_",
                                    var_idx
                            );
                        }
                    }

                    int l_reach = 0, r_reach = 0, score = 0;
                    if (left_compute || right_compute) {
                        score = calc_vcf_swg_score(
                                vars, clust, clust+1, sub, open, extend);
                        if (print) printf("    orig score: %d\n", score);
                    }

                    // LEFT REACH 

                    if (left_compute) { // calculate left reach
                        std::string query, ref;

                        // error checking
                        if (vars->clusters[clust]-1 < 0)
                            ERROR("left var_idx < 0");
                        if (vars->clusters[clust]-1 >= nvar)
                            ERROR("left var_idx >= nvar");
                        if (clust+1 >= vars->clusters.size())
                            ERROR("left next clust_idx >= nclust");
                        if (vars->clusters[clust+1]-1 < 0)
                            ERROR("left next var_idx < 0");
                        if (vars->clusters[clust+1]-1 >= nvar)
                            ERROR("left next var_idx >= nvar");

                        // just after last variant in this cluster
                        int end_pos = vars->poss[vars->clusters[clust+1]-1] +
                                vars->rlens[vars->clusters[clust+1]-1]+1;
                        
                        // get reference end pos of last variant in this cluster
                        int main_diag_start = end_pos - vars->poss[ 
                                    vars->clusters[clust]];
                        int main_diag = 0;
                        for (int vi = vars->clusters[clust];
                                vi < vars->clusters[clust+1]; vi++)
                            main_diag += vars->refs[vi].size() - vars->alts[vi].size();

                        // calculate max reaching path to left
                        int ref_len = score/extend + 3;
                        int reach = ref_len - 1;
                        while (reach == ref_len-1) {
                            ref_len *= 2;
                            int beg_pos = end_pos - std::max(0, main_diag) - ref_len - score/extend-3;
g.timers[TIME_CL_GENSTR].start();
                            query = generate_str(vcf->ref, vars, ctg, 
                                        vars->clusters[clust], 
                                        vars->clusters[clust+1], 
                                        beg_pos, end_pos);
g.timers[TIME_CL_GENSTR].stop();
g.timers[TIME_CL_SUBSTR].start();
                            ref = vcf->ref->fasta.at(ctg).substr(end_pos-ref_len, ref_len);
                            std::reverse(query.begin(), query.end());
                            std::reverse(ref.begin(), ref.end());
g.timers[TIME_CL_SUBSTR].stop();
g.timers[TIME_CL_REACH].start();
g.timers[TIME_CL_BUFF_INIT].start();
                            // manage buffer for storing offsets
                            size_t offs_size = MATS * (std::max(open+extend, sub)+1) * 
                                (query.size() + ref.size() - 1);
                            if (offs_size > offs_buffer.size())
                                offs_buffer.resize(offs_size, -2);
g.timers[TIME_CL_BUFF_INIT].stop();
                            // calculate reach
                            reach = wf_swg_max_reach(query, ref, offs_buffer,
                                    main_diag, main_diag_start, score,
                                    sub, open, extend, 
                                    false /* print */, true  /* reverse */);
                            // reset buffer
g.timers[TIME_CL_BUFF_CLEAR].start();
                            for (size_t i = 0; i < offs_size; i++)
                                offs_buffer[i] = -2;
g.timers[TIME_CL_BUFF_CLEAR].stop();
g.timers[TIME_CL_REACH].stop();
                        }
                        if (print) printf("    left reach: %d\n", reach);
                        l_reach = end_pos - reach;

                        if (false) {
                            printf("REF:        %s\n", ref.data());
                            printf("QUERY:      %s\n", query.data());
                            printf(" main diag:  %d\n", main_diag);
                            printf("diag start:  %d\n\n", main_diag_start);
                        }
                        
                    } else {
                        // past farthest right (unused)
                        if (clust < prev_clusters.size()-1)
                            l_reach = vars->poss[vars->clusters[clust]];
                        else
                            l_reach = len + g.reach_min_gap*2;
                    }
                    left_reach.push_back(l_reach);

                    // RIGHT REACH
                    if (right_compute) { // calculate right reach
                        std::string query, ref;

                        // error checking
                        if (vars->clusters[clust] < 0)
                            ERROR("right var_idx < 0");
                        if (vars->clusters[clust] >= nvar)
                            ERROR("right var_idx >= nvar");
                        if (clust+1 >= vars->clusters.size())
                            ERROR("right next clust_idx >= nclust");
                        if (vars->clusters[clust+1] < 0)
                            ERROR("right next var_idx < 0");
                        if (vars->clusters[clust+1] >= nvar)
                            ERROR("right next var_idx >= nvar");

                        // right before current cluster
                        int beg_pos = vars->poss[vars->clusters[clust]]-1;

                        // get reference end pos of last variant in this cluster
                        int main_diag_start = vars->poss[ 
                                    vars->clusters[clust+1]-1]
                                    + vars->rlens[vars->clusters[clust+1]-1] - beg_pos;
                        int main_diag = 0;
                        for (int vi = vars->clusters[clust];
                                vi < vars->clusters[clust+1]; vi++)
                            main_diag += vars->refs[vi].size() - vars->alts[vi].size();

                        // calculate max reaching path to right
                        int ref_len = score/extend + 3;
                        int reach = ref_len - 1;
                        while (reach == ref_len-1) {
                            ref_len *= 2;
                            int end_pos = beg_pos + std::max(0, main_diag) + ref_len + score/extend+3;
g.timers[TIME_CL_GENSTR].start();
                            query = generate_str(vcf->ref, vars, ctg, 
                                        vars->clusters[clust], 
                                        vars->clusters[clust+1], 
                                        beg_pos, end_pos);
g.timers[TIME_CL_GENSTR].stop();
g.timers[TIME_CL_SUBSTR].start();
                            ref = vcf->ref->fasta.at(ctg).substr(beg_pos, ref_len);
g.timers[TIME_CL_SUBSTR].stop();
g.timers[TIME_CL_REACH].start();
                            // manage buffer for storing offsets
                            size_t offs_size = MATS * (std::max(sub, open+extend)+1) * 
                                (query.size() + ref.size() - 1);
                            if (offs_size > offs_buffer.size())
                                offs_buffer.resize(offs_size, -2);
                            // calculate reach
                            reach = wf_swg_max_reach(query, ref, offs_buffer,
                                    main_diag, main_diag_start, score,
                                    sub, open, extend, 
                                    false /* print */, false /* reverse */);
                            // reset buffer
                            for (size_t i = 0; i < offs_size; i++)
                                offs_buffer[i] = -2;
g.timers[TIME_CL_REACH].stop();
                        }
                        if (print) printf("   right reach: %d\n", reach);
                        r_reach = beg_pos + reach + 1;

                        if (false) {
                            printf("REF:        %s\n", ref.data());
                            printf("QUERY:      %s\n", query.data());
                            printf(" main diag:  %d\n", main_diag);
                            printf("diag start:  %d\n\n", main_diag_start);
                        }
                        
                    } else { // non-adjacent, don't really compute
                        if (clust < prev_clusters.size()-1)
                            r_reach = vars->poss[vars->clusters[clust+1]-1] +
                                     vars->rlens[vars->clusters[clust+1]-1];
                        else
                            r_reach = -g.reach_min_gap*2; // past farthest left (unused)
                    }
                    right_reach.push_back(r_reach);
                    if (print) printf("span: %s - %s\n", 
                                l_reach == len+g.reach_min_gap*2 ? 
                                    "X" : std::to_string(l_reach).data(), 
                                r_reach == -g.reach_min_gap*2 ? "X" : std::to_string(r_reach).data());
                }

                // merge dependent clusters rightwards
                std::vector<int> tmp_left_reach, tmp_right_reach;
                tmp_right_reach.push_back(-g.reach_min_gap*2);
                tmp_left_reach.push_back(len+g.reach_min_gap*2);
                int clust = 0;
                while (clust < int(prev_clusters.size())) {
                    int clust_size = 1;
                    int max_right_reach = right_reach[clust];
                    int min_left_reach = left_reach[clust];
                    while (clust+clust_size < int(prev_clusters.size()) &&
                            max_right_reach + g.reach_min_gap >= left_reach[clust+clust_size]) {
                        // keep comparing against farthest-right seen so far
                        max_right_reach = std::max(max_right_reach,
                                right_reach[clust+clust_size]);
                        // calculate to save for next step
                        min_left_reach = std::min(min_left_reach,
                                left_reach[clust+clust_size]);
                        clust_size++;
                    }
                    tmp_right_reach.push_back(max_right_reach);
                    tmp_left_reach.push_back(min_left_reach);
                    tmp_clusters.push_back(prev_clusters[clust]);
                    tmp_active.push_back(clust_size > 1); // true if merge occurred
                    clust += clust_size;
                }

                /* printf("tmp clusters:"); */
                /* for (int i = 0; i < int(tmp_clusters.size()); i++) { */
                /*     printf("\t%d", tmp_clusters[i]); */
                /* } */
                /* printf("\n"); */
                /* printf("  tmp active:"); */
                /* for (int i = 0; i < int(tmp_active.size()); i++) { */
                /*     printf("\t%d", tmp_active[i] ? 1 : 0); */
                /* } */
                /* printf("\n"); */

                // merge dependent clusters leftwards
                clust = tmp_clusters.size()-1;
                bool was_merged = false;
                while (clust >= 0) {
                    int clust_size = 1;
                    int min_left_reach = tmp_left_reach[clust];
                    while (clust - clust_size >= 0 &&
                            min_left_reach <= 
                            tmp_right_reach[clust-clust_size] + g.reach_min_gap) {
                        min_left_reach = std::min(min_left_reach,
                                tmp_left_reach[clust-clust_size]);
                        clust_size++;
                    }
                    next_clusters.push_back(tmp_clusters[clust]);
                    next_active.push_back(was_merged); // from previous itr
                    was_merged = clust_size > 1 || tmp_active[clust];
                    clust -= clust_size;
                }
                std::reverse(next_clusters.begin(), next_clusters.end());
                std::reverse(next_active.begin(), next_active.end());

                /* printf("next clusters:"); */
                /* for (int i = 0; i < int(next_clusters.size()); i++) { */
                /*     printf("\t%d", next_clusters[i]); */
                /* } */
                /* printf("\n"); */
                /* printf("  next active:"); */
                /* for (int i = 0; i < int(next_active.size()); i++) { */
                /*     printf("\t%d", next_active[i] ? 1 : 0); */
                /* } */
                /* printf("\n"); */

                // reset for next iteration
                prev_clusters = next_clusters;
                prev_active = next_active;
                tmp_clusters.clear();
                tmp_active.clear();
                next_clusters.clear();
                next_active.clear();
                left_reach.clear();
                right_reach.clear();
            }

            // save final clustering
            vars->clusters = prev_clusters;

            if (nvar && g.verbosity >= 1) INFO("      %d resulting clusters.", int(prev_clusters.size()-1));
        }
    }
}
