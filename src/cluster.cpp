#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "globals.h"
#include "dist.h"
#include "cluster.h"
#include "print.h"
#include "variant.h"


ctgSuperclusters::ctgSuperclusters() {
    this->callset_vars = std::vector< std::shared_ptr<ctgVariants> > (CALLSETS, nullptr);
}


/******************************************************************************/


int ctgSuperclusters::get_min_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end) {
    // NOTE: qvi_start == qvi_end == vars->n is valid (empty sentinel), should return int::max
    return std::min( (qvi_start == qvi_end) ? std::numeric_limits<int>::max() : 
                this->callset_vars[QUERY]->poss[qvi_start], 
            (tvi_start == tvi_end) ? std::numeric_limits<int>::max() :
                this->callset_vars[TRUTH]->poss[tvi_start]) - 1;
}


int ctgSuperclusters::get_max_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end) {
    // NOTE: qvi_start == qvi_end == vars->n is valid (empty sentinel), should return int::max
    int max_ref_pos = -1;
    for (int qvi = qvi_start; qvi < qvi_end; qvi++) {
        max_ref_pos = std::max(max_ref_pos,
                this->callset_vars[QUERY]->poss[qvi] + this->callset_vars[QUERY]->rlens[qvi]);
    }
    for (int tvi = tvi_start; tvi < tvi_end; tvi++) {
        max_ref_pos = std::max(max_ref_pos,
                this->callset_vars[TRUTH]->poss[tvi] + this->callset_vars[TRUTH]->rlens[tvi]);
    }
    return (max_ref_pos == -1) ? std::numeric_limits<int>::max() : max_ref_pos + 1;
}


/******************************************************************************/


std::vector< std::vector< std::vector<int> > > 
sort_superclusters(std::shared_ptr<superclusterData> sc_data) {

    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Sorting superclusters by size");
    std::vector< std::vector< std::vector<int> > > sc_groups(g.thread_nsteps,
            std::vector< std::vector<int> >(2));

    for (int ctg_idx = 0; ctg_idx < int(sc_data->contigs.size()); ctg_idx++) {
        std::string ctg = sc_data->contigs[ctg_idx];
        std::shared_ptr<ctgSuperclusters> ctg_scs = sc_data->superclusters[ctg];
        std::shared_ptr<ctgVariants> qvars = sc_data->superclusters[ctg]->callset_vars[QUERY];
        std::shared_ptr<ctgVariants> tvars = sc_data->superclusters[ctg]->callset_vars[TRUTH];
        if (!qvars->n) continue;
        // superclusters are numbered 0...n-1, so we need +1 to get the total count
        int nscs =  std::max(qvars->superclusters[qvars->n-1], tvars->superclusters[tvars->n-1]) + 1;

        for (int sc_idx = 0; sc_idx < nscs; sc_idx++) {

            std::vector<size_t> max_lens(CALLSETS, 0);
            std::vector<size_t> lens(HAPS, 0);
            for (int c = 0; c < CALLSETS; c++) {
                auto vars = ctg_scs->callset_vars[c];
                if (ctg_scs->callset_vars[c]->nc == 0) continue;
                int var_beg = std::distance(vars->superclusters.begin(),
                    std::lower_bound(vars->superclusters.begin(), vars->superclusters.end(), sc_idx));
                int var_end = std::distance(vars->superclusters.begin(),
                    std::upper_bound(vars->superclusters.begin(), vars->superclusters.end(), sc_idx));

                // calculate query len as reference length plus alternate lengths
                for (int hi = 0; hi < HAPS; hi++) {
                    lens[hi] = var_beg == var_end ? 0 : vars->poss[var_end-1] - vars->poss[var_beg];
                    for (int vi = var_beg; vi < var_end; vi++) {
                        if (!vars->var_on_hap(vi, hi)) continue;
                        lens[hi] += int(vars->alts[vi].size());
                    }
                    max_lens[c] = std::max(max_lens[c], lens[hi]);
                }
            }

            // calculate memory usage
            int FRAG_FACTOR = 2;
            size_t mem = (max_lens[QUERY] + max_lens[TRUTH]) * (sizeof(uint32_t) + sizeof(int)) * 
                g.max_dist * FRAG_FACTOR;
            double mem_gb = mem / (1000.0 * 1000.0 * 1000.0);
            if (mem_gb > g.max_ram) {
                WARN("Max (%.3fGB) RAM exceeded (%.3fGB req) for supercluster %d, running anyways. Lower --max-supercluster-size if this fails.", 
                        g.max_ram, mem_gb, sc_idx);
                sc_groups[g.thread_nsteps-1][CTG_IDX].push_back(ctg_idx);
                sc_groups[g.thread_nsteps-1][SC_IDX].push_back(sc_idx);
                continue;
            }
            
            // place into correct group
            for (int i = 0; i < g.thread_nsteps; i++) {
                if (mem_gb < g.ram_steps[i]) {
                    sc_groups[i][CTG_IDX].push_back(ctg_idx);
                    sc_groups[i][SC_IDX].push_back(sc_idx);
                    break;
                }
            }
        }
    }

    return sc_groups;
}

/******************************************************************************/

void superclusterData::add_callset_vars(int callset,
        std::vector< std::unordered_map< std::string, 
        std::shared_ptr<ctgVariants> > > & vars) {
    bool print = false;

    for (int ctg_idx = 0; ctg_idx < int(this->contigs.size()); ctg_idx++) {
        std::string ctg = contigs[ctg_idx];

        // merge variants from each haplotype for this ctg/callset
        std::shared_ptr<ctgVariants> merged_vars(new ctgVariants(ctg));

        // skip empty contigs
        int nvars = 0;
        for (int h = 0; h < HAPS; h++) nvars += vars[h][ctg]->n;
        if (!nvars) {
            this->superclusters[ctg]->callset_vars[callset] = merged_vars;
            continue;
        }

        // initialize with leftmost cluster info
        int curr_left_reach = std::numeric_limits<int>::max();
        int curr_right_reach = std::numeric_limits<int>::min();
        if (vars[HAP1][ctg]->n && vars[HAP2][ctg]->n) {
            if (vars[HAP1][ctg]->left_reaches[0] < vars[HAP2][ctg]->left_reaches[0]) {
                curr_left_reach = vars[HAP1][ctg]->left_reaches[0];
                curr_right_reach = vars[HAP1][ctg]->right_reaches[0];
            } else {
                curr_left_reach = vars[HAP2][ctg]->left_reaches[0];
                curr_right_reach = vars[HAP2][ctg]->right_reaches[0];
            }
        } else if (vars[HAP1][ctg]->n) {
            curr_left_reach = vars[HAP1][ctg]->left_reaches[0];
            curr_right_reach = vars[HAP1][ctg]->right_reaches[0];
        } else if (vars[HAP2][ctg]->n) {
            curr_left_reach = vars[HAP2][ctg]->left_reaches[0];
            curr_right_reach = vars[HAP2][ctg]->right_reaches[0];
        } else {
            ERROR("No variants on contig: %s", ctg.data());
        }

        // initialize indices
        std::vector<int> var_idx(HAPS, 0);
        std::vector<int> clust_idx(HAPS, 0);
        int curr_var_idx = 0;
        int next_left_reach = std::numeric_limits<int>::max();
        int next_right_reach = std::numeric_limits<int>::min();
        while (var_idx[HAP1] < vars[HAP1][ctg]->n || var_idx[HAP2] < vars[HAP2][ctg]->n) {

            //////////////////
            // ADD VARIANTS //
            //////////////////
            std::vector<bool> updated(HAPS, false);

            if (var_idx[HAP1] < vars[HAP1][ctg]->n &&
                    var_idx[HAP2] < vars[HAP2][ctg]->n) {

                // variants are same, add homozygous
                if (vars[HAP1][ctg]->poss[var_idx[HAP1]] == 
                        vars[HAP2][ctg]->poss[var_idx[HAP2]] && 
                        vars[HAP1][ctg]->refs[var_idx[HAP1]] ==
                        vars[HAP2][ctg]->refs[var_idx[HAP2]] &&
                        vars[HAP1][ctg]->alts[var_idx[HAP1]] ==
                        vars[HAP2][ctg]->alts[var_idx[HAP2]]) {
                    merged_vars->add_var(
                            vars[HAP1][ctg]->poss[var_idx[HAP1]],
                            vars[HAP1][ctg]->rlens[var_idx[HAP1]],
                            vars[HAP1][ctg]->types[var_idx[HAP1]],
                            vars[HAP1][ctg]->locs[var_idx[HAP1]],
                            vars[HAP1][ctg]->refs[var_idx[HAP1]],
                            vars[HAP1][ctg]->alts[var_idx[HAP1]],
                            GT_ALT1_ALT1,
                            vars[HAP1][ctg]->gt_quals[var_idx[HAP1]],
                            vars[HAP1][ctg]->var_quals[var_idx[HAP1]],
                            vars[HAP1][ctg]->phase_sets[var_idx[HAP1]]);
                    if (print) printf("adding 1|1 var= %s:%d\t%s\t%s\t%s\n",
                            ctg.data(), 
                            vars[HAP1][ctg]->poss[var_idx[HAP1]],
                            vars[HAP1][ctg]->refs[var_idx[HAP1]].data(),
                            vars[HAP1][ctg]->alts[var_idx[HAP1]].data(),
                            gt_strs[vars[HAP1][ctg]->orig_gts[var_idx[HAP1]]].data()
                            );
                    var_idx[HAP1]++; var_idx[HAP2]++;
                    updated[HAP1] = true; updated[HAP2] = true;

                } else { // heterozygous, add first-occurring variant
                    int hap_idx = HAP1; // default to HAP1
                    // if location is same, default to INS first
                    if (vars[HAP1][ctg]->poss[var_idx[HAP1]] ==
                            vars[HAP2][ctg]->poss[var_idx[HAP2]] &&
                            vars[HAP2][ctg]->types[var_idx[HAP2]] == TYPE_INS) {
                        hap_idx = HAP2;
                    // otherwise, choose first variant
                    } else if (vars[HAP2][ctg]->poss[var_idx[HAP2]] <
                            vars[HAP1][ctg]->poss[var_idx[HAP1]]) {
                        hap_idx = HAP2;
                    }
                    merged_vars->add_var(
                            vars[hap_idx][ctg]->poss[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->rlens[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->types[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->locs[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->refs[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->alts[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->orig_gts[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->gt_quals[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->var_quals[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->phase_sets[var_idx[hap_idx]]);
                    if (print) printf("adding %s var= %s:%d\t%s\t%s\t%s\n",
                            hap_idx == 0 ? "1|0" : "0|1",
                            ctg.data(), 
                            vars[hap_idx][ctg]->poss[var_idx[hap_idx]],
                            vars[hap_idx][ctg]->refs[var_idx[hap_idx]].data(),
                            vars[hap_idx][ctg]->alts[var_idx[hap_idx]].data(),
                            gt_strs[vars[hap_idx][ctg]->orig_gts[var_idx[hap_idx]]].data()
                                );
                    var_idx[hap_idx]++;
                    updated[hap_idx] = true;
                }
            } else if (var_idx[HAP1] < vars[HAP1][ctg]->n) { // only hap1 vars left
                merged_vars->add_var(
                        vars[HAP1][ctg]->poss[var_idx[HAP1]],
                        vars[HAP1][ctg]->rlens[var_idx[HAP1]],
                        vars[HAP1][ctg]->types[var_idx[HAP1]],
                        vars[HAP1][ctg]->locs[var_idx[HAP1]],
                        vars[HAP1][ctg]->refs[var_idx[HAP1]],
                        vars[HAP1][ctg]->alts[var_idx[HAP1]],
                        vars[HAP1][ctg]->orig_gts[var_idx[HAP1]],
                        vars[HAP1][ctg]->gt_quals[var_idx[HAP1]],
                        vars[HAP1][ctg]->var_quals[var_idx[HAP1]],
                        vars[HAP1][ctg]->phase_sets[var_idx[HAP1]]);
                if (print) printf("adding 1|0 var= %s:%d\t%s\t%s\t%s\n",
                        ctg.data(), 
                        vars[HAP1][ctg]->poss[var_idx[HAP1]],
                        vars[HAP1][ctg]->refs[var_idx[HAP1]].data(),
                        vars[HAP1][ctg]->alts[var_idx[HAP1]].data(),
                        gt_strs[vars[HAP1][ctg]->orig_gts[var_idx[HAP1]]].data()
                            );
                var_idx[HAP1]++;
                updated[HAP1] = true;
            } else if (var_idx[HAP2] < vars[HAP2][ctg]->n) { // only hap2 vars left
                merged_vars->add_var(
                        vars[HAP2][ctg]->poss[var_idx[HAP2]],
                        vars[HAP2][ctg]->rlens[var_idx[HAP2]],
                        vars[HAP2][ctg]->types[var_idx[HAP2]],
                        vars[HAP2][ctg]->locs[var_idx[HAP2]],
                        vars[HAP2][ctg]->refs[var_idx[HAP2]],
                        vars[HAP2][ctg]->alts[var_idx[HAP2]],
                        vars[HAP2][ctg]->orig_gts[var_idx[HAP2]],
                        vars[HAP2][ctg]->gt_quals[var_idx[HAP2]],
                        vars[HAP2][ctg]->var_quals[var_idx[HAP2]],
                        vars[HAP2][ctg]->phase_sets[var_idx[HAP2]]);
                if (print) printf("adding 0|1 var= %s:%d\t%s\t%s\t%s\n",
                        ctg.data(), 
                        vars[HAP2][ctg]->poss[var_idx[HAP2]],
                        vars[HAP2][ctg]->refs[var_idx[HAP2]].data(),
                        vars[HAP2][ctg]->alts[var_idx[HAP2]].data(),
                        gt_strs[vars[HAP2][ctg]->orig_gts[var_idx[HAP2]]].data());
                var_idx[HAP2]++;
                updated[HAP2] = true;
            }

            /////////////////////
            // UPDATE CLUSTERS //
            /////////////////////

            // update cluster indices before checking reaches
            for (int h = 0; h < HAPS; h++) {
                if (updated[h] && clust_idx[h] < vars[h][ctg]->nc) {
                    if (var_idx[h] >= vars[h][ctg]->clusters[clust_idx[h]+1]) {
                        clust_idx[h]++;
                    }
                }
            }

            // update reaches
            for (int h = 0; h < HAPS; h++) {
                if (clust_idx[h] < vars[h][ctg]->nc) {
                    if (clust_idx[h^1] >= vars[h^1][ctg]->nc ||
                            vars[h][ctg]->left_reaches[clust_idx[h]] <=
                            vars[h^1][ctg]->left_reaches[clust_idx[h^1]]) {
                        // TODO: why no std::min here?
                        next_left_reach = vars[h][ctg]->left_reaches[clust_idx[h]];
                        next_right_reach = std::max(next_right_reach,
                            vars[h][ctg]->right_reaches[clust_idx[h]]);
                    }
                }
            }
            if (print) printf("curr = (%d, %d)\tnext = (%d, %d)\n",
                    curr_left_reach, curr_right_reach,
                    next_left_reach, next_right_reach);

            // starting new supercluster
            if (next_left_reach > curr_right_reach) {
                if (print) printf("new supercluster, adding curr = %d (%d, %d)\n",
                        curr_var_idx, curr_left_reach, curr_right_reach);

                // save reaches of prev cluster
                merged_vars->clusters.push_back(curr_var_idx);
                merged_vars->left_reaches.push_back(curr_left_reach);
                merged_vars->right_reaches.push_back(curr_right_reach);
                merged_vars->nc++;

                // set reaches of this cluster
                curr_var_idx = merged_vars->n;
                curr_left_reach  = next_left_reach;
                curr_right_reach = next_right_reach;

                // init reaches of next cluster
                next_left_reach  = std::numeric_limits<int>::max();
                next_right_reach = std::numeric_limits<int>::min();
                for (int h = 0; h < HAPS; h++) {
                    if (clust_idx[h] < vars[h][ctg]->nc &&
                       (clust_idx[h^1] >= vars[h^1][ctg]->nc ||
                            vars[h][ctg]->left_reaches[clust_idx[h]] <=
                            vars[h^1][ctg]->left_reaches[clust_idx[h^1]])) {
                        next_left_reach  = std::min(next_left_reach,
                                vars[h][ctg]->left_reaches[clust_idx[h]]);
                        next_right_reach = std::max(next_right_reach,
                                vars[h][ctg]->right_reaches[clust_idx[h]]);
                    }
                }

            } else { // same supercluster
                if (print) printf("same supercluster\n");
                curr_left_reach = std::min(curr_left_reach, next_left_reach);
                curr_right_reach = std::max(curr_right_reach, next_right_reach);
            }

        } // while variants remain

        // save reaches of final cluster
        merged_vars->clusters.push_back(curr_var_idx);
        merged_vars->left_reaches.push_back(curr_left_reach);
        merged_vars->right_reaches.push_back(curr_right_reach);
        merged_vars->nc++;

        // add sentinel cluster and save contig variants
        merged_vars->clusters.push_back(merged_vars->n);
        merged_vars->left_reaches.push_back(std::numeric_limits<int>::max());
        merged_vars->right_reaches.push_back(std::numeric_limits<int>::max());
        merged_vars->nc++;
        this->superclusters[ctg]->callset_vars[callset] = merged_vars;

    } // for each contig
}

/******************************************************************************/

superclusterData::superclusterData(
        std::shared_ptr<variantData> query_ptr,
        std::shared_ptr<variantData> truth_ptr,
        std::shared_ptr<fastaData> ref_ptr) {

    // set reference pointer
    this->ref = ref_ptr;

    // save samples and filenames
    this->samples.push_back(query_ptr->sample);
    this->samples.push_back(truth_ptr->sample);
    this->filenames.push_back(query_ptr->filename);
    this->filenames.push_back(truth_ptr->filename);

    // create list of all contigs covered by truth/query
    for (int i = 0; i < int(query_ptr->contigs.size()); i++) {
        std::string ctg = query_ptr->contigs[i];
        if (std::find(this->contigs.begin(), this->contigs.end(), ctg) == 
                this->contigs.end()) {
            this->contigs.push_back(ctg);
            this->lengths.push_back(query_ptr->lengths[i]);
            this->ploidy.push_back(query_ptr->ploidy[i]);
            this->superclusters[ctg] = std::shared_ptr<ctgSuperclusters>(
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
            this->superclusters[ctg] = std::shared_ptr<ctgSuperclusters>(
                    new ctgSuperclusters());
        }
   }

    // for each contig, merge variants and clusters across haplotypes
    this->add_callset_vars(QUERY, query_ptr->variants);
    this->add_callset_vars(TRUTH, truth_ptr->variants);
}


/******************************************************************************/

/* Supercluster using left_reach and right_reach of each cluster, calculated
 * during wf_swg_cluster().
 */
void superclusterData::supercluster(bool print) {

    // iterate over each contig
    int total_superclusters = 0;
    int largest_supercluster = 0;
    int total_vars = 0;
    int most_vars = 0;
    int total_bases = 0;
    for (const std::string & ctg : this->contigs) {

        // skip empty contigs
        int nvars = 0;
        for (int c = 0; c < CALLSETS; c++) {
            nvars += this->superclusters[ctg]->callset_vars[c]->n;
        }
        if (!nvars) continue;

        // for each cluster of variants (merge query and truth haps)
        auto & vars = this->superclusters[ctg]->callset_vars;
        std::vector<int> brks(CALLSETS, 0); // start of current supercluster
        int sc_idx = 0;
        while (true) {

            // init: empty supercluster
            std::vector<int> next_brks = brks; // end of current supercluster
            std::vector<int> lefts(CALLSETS, std::numeric_limits<int>::max());
            for (int c = 0; c < CALLSETS; c++) {
                if (brks[c] < vars[c]->nc) {
                    lefts[c] = vars[c]->left_reaches[ next_brks[c] ];
                }
            }

            // get first cluster, end if all haps are off end
            int c = std::distance(lefts.begin(),
                    std::min_element(lefts.begin(), lefts.end()));
            if (lefts[c] == std::numeric_limits<int>::max()) break;

            // initialize cluster merging with first to start
            int curr_right = vars[c]->right_reaches[ next_brks[c] ];
            next_brks[c]++;
            lefts[c] = next_brks[c] < vars[c]->nc ?
                    vars[c]->left_reaches[ next_brks[c] ] :
                    std::numeric_limits<int>::max();

            // keep expanding cluster while possible
            bool just_active = true;
            while (just_active) {
                just_active = false;
                for (int c = 0; c < CALLSETS; c++) {
                    while (lefts[c] <= curr_right) {
                        curr_right = std::max(curr_right,
                            vars[c]->right_reaches[ next_brks[c] ]);
                        next_brks[c]++;
                        lefts[c] = next_brks[c] < vars[c]->nc ?
                            vars[c]->left_reaches[ next_brks[c] ] :
                            std::numeric_limits<int>::max();
                        just_active = true;
                    }
                }
            }

            // get the range of positions spanned by this supercluster
            std::vector<int> poss = get_supercluster_range(vars, brks, next_brks);
            int beg_pos = poss[0];
            int end_pos = poss[1];

            // debug print
            if (print) {
                printf("\nSUPERCLUSTER: %d\n", sc_idx);
                printf("POS: %s:%d-%d\n", ctg.data(), beg_pos, end_pos);
                printf("SIZE: %d\n", end_pos - beg_pos);
                for (int ci = 0; ci < CALLSETS; ci++) {
                    printf("%s: clusters %d-%d of %d = %d\n",
                            callset_strs[ci].data(),
                            brks[ci], next_brks[ci], vars[ci]->nc, int(vars[ci]->clusters.size()));
                    if (vars[ci]->nc)
                    for (int vi = vars[ci]->clusters[ brks[ci]]; 
                            vi < vars[ci]->clusters[ next_brks[ci] ]; vi++) {
                        printf("    var %d of %d = %s:%d\t%s\t%s\t%s\n", vi, vars[ci]->n, ctg.data(),
                                vars[ci]->poss[vi], 
                                vars[ci]->refs[vi].data(),
                                vars[ci]->alts[vi].data(),
                                gt_strs[vars[ci]->orig_gts[vi]].data());
                    }
                }
                printf("curr_right: %d, lefts: %d %d\n", curr_right, lefts[0], lefts[1]);
            }

            // split large supercluster if necessary
            if (end_pos - beg_pos > g.max_supercluster_size) {
                std::vector< std::vector<int> > all_brks = split_large_supercluster(vars, brks, next_brks, print);
                WARN("Max supercluster size (%d) exceeded (%d) at %s:%d-%d, breaking up into %d superclusters",
                        g.max_supercluster_size, end_pos - beg_pos, ctg.data(), beg_pos, end_pos, int(all_brks.size()-1));
                for (int brk_idx = 0; brk_idx < int(all_brks.size())-1; brk_idx++) {
                    brks = all_brks[brk_idx];
                    next_brks = all_brks[brk_idx+1];
                    poss = get_supercluster_range(vars, brks, next_brks);
                    beg_pos = poss[0];
                    end_pos = poss[1];

                    // save supercluster information
                    for (int ci = 0; ci < CALLSETS; ci++) {
                        for (int vi = vars[ci]->clusters[brks[ci]]; 
                                 vi < vars[ci]->clusters[next_brks[ci]]; vi++) {
                            vars[ci]->superclusters[vi] = sc_idx;
                        }
                    }
                    sc_idx++;

                    // update summary metrics
                    largest_supercluster = std::max(largest_supercluster, end_pos-beg_pos);
                    total_bases += end_pos-beg_pos;
                    int this_vars = 0;
                    for (int ci = 0; ci < CALLSETS; ci++) {
                        if (vars[ci]->nc)
                            this_vars += vars[ci]->clusters[next_brks[ci]] - vars[ci]->clusters[brks[ci]];
                    }
                    most_vars = std::max(most_vars, this_vars);
                    total_vars += this_vars;
                }
            } else {

                // save supercluster information
                for (int ci = 0; ci < CALLSETS; ci++) {
                    for (int vi = vars[ci]->clusters[brks[ci]]; 
                             vi < vars[ci]->clusters[next_brks[ci]]; vi++) {
                        vars[ci]->superclusters[vi] = sc_idx;
                    }
                }
                sc_idx++;

                // update summary metrics
                largest_supercluster = std::max(largest_supercluster, end_pos-beg_pos);
                total_bases += end_pos-beg_pos;
                int this_vars = 0;
                for (int ci = 0; ci < CALLSETS; ci++) {
                    if (vars[ci]->nc)
                        this_vars += vars[ci]->clusters[next_brks[ci]] - vars[ci]->clusters[brks[ci]];
                }
                most_vars = std::max(most_vars, this_vars);
                total_vars += this_vars;

            }

            // reset for next active cluster
            brks = next_brks;
        }

        // add remaining variants as a supercluster
        for (int ci = 0; ci < CALLSETS; ci++) {
            for (int vi = vars[ci]->clusters[brks[ci]]; vi < vars[ci]->n; vi++) {
                vars[ci]->superclusters[vi] = sc_idx;
            }
        }

        total_superclusters += sc_idx;
    }
    if (g.verbosity >= 1) INFO("           Total superclusters: %d", total_superclusters);
    if (g.verbosity >= 1) INFO("  Largest supercluster (bases): %d", largest_supercluster);
    if (g.verbosity >= 1) INFO("  Largest supercluster  (vars): %d", most_vars);
    if (g.verbosity >= 1 && total_superclusters) INFO("  Average supercluster (bases): %.3f", total_bases / float(total_superclusters));
    if (g.verbosity >= 1 && total_superclusters) INFO("  Average supercluster  (vars): %.3f", total_vars / float(total_superclusters));
}


/******************************************************************************/


/* Calculate the size of a supercluster from a list of variants and the start/stop indices of the
 * clusters that compose the supercluster. Returns a vector containing the (start, end) positions.
 */
std::vector<int> get_supercluster_range(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices, // inclusive
        const std::vector<int> & cluster_end_indices) { // exclusive

    int beg_pos = std::numeric_limits<int>::max();
    int end_pos = -1;

    for (int c = 0; c < CALLSETS; c++) {
        // if there is a cluster on this hap, update beginning and end positions
        if (cluster_end_indices[c] - cluster_start_indices[c]) {
            if (cluster_start_indices[c] > vars[c]->nc) {
                ERROR("Cluster start indices invalid in get_supercluster_range()");
            }
            if (cluster_end_indices[c] > vars[c]->nc) {
                ERROR("Cluster end indices invalid in get_supercluster_range()");
            }

            // one position left of the leftmost variant
            beg_pos = std::min(beg_pos,
                vars[c]->poss[vars[c]->clusters[cluster_start_indices[c]]]-1);

            // one position right of the rightmost variant end in the last included cluster
            end_pos = std::max(end_pos,
                    vars[c]->poss[vars[c]->clusters[cluster_end_indices[c]]-1] +
                    vars[c]->rlens[vars[c]->clusters[cluster_end_indices[c]]-1] + 1);
        }
    }
    std::vector<int> range = {beg_pos, end_pos};
    return range;
}


/******************************************************************************/


/* Split the supercluster as many times as necessary until it is within the size limits.
 * Note: this function is only called when initial supercluster is too large.
 */
std::vector< std::vector<int> > split_large_supercluster(
        std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        std::vector<int> & cluster_end_indices, bool print) {

    std::vector< std::vector<int> > breakpoints = {cluster_start_indices, cluster_end_indices};
    bool large_supercluster_exists = true;
    while (large_supercluster_exists) {
        large_supercluster_exists = false;
        
        std::vector< std::vector<int> > next_breakpoints;
        for (int i = 0; i < int(breakpoints.size())-1; i++) {
            std::vector<int> poss = get_supercluster_range(vars, breakpoints[i], breakpoints[i+1]);
            int beg_pos = poss[0];
            int end_pos = poss[1];
            if (end_pos - beg_pos > g.max_supercluster_size) {
                if (print) printf("\nsplitting %d-%d\n", beg_pos, end_pos);
                large_supercluster_exists = true;
                next_breakpoints.push_back(breakpoints[i]);

                // variant indices of optimal split location
                std::vector<int> best_var_split = get_supercluster_split_location(
                            vars, breakpoints[i], breakpoints[i+1], print);

                if (int(best_var_split.size()) == CALLSETS) { // found valid split
                    std::vector<int> cluster_split_indices = split_cluster(vars, best_var_split, breakpoints, i, print);
                    next_breakpoints.push_back(cluster_split_indices);

                } else { // no valid splits (shouldn't happen?)
                    printf("WARNING: no valid splits at %d-%d\n", beg_pos, end_pos);
                    large_supercluster_exists = false;
                }

            } else { // this supercluster is small enough
                next_breakpoints.push_back(breakpoints[i]);
            }
        }
        // add last breakpoint
        cluster_end_indices = breakpoints[breakpoints.size()-1];
        next_breakpoints.push_back(cluster_end_indices);
        breakpoints = next_breakpoints;
    }
    if (print) printf("%d breakpoints\n", int(breakpoints.size()));
     
    return breakpoints;
}


/******************************************************************************/

/* Given a 2-tuple of variant indices (on each callset), at which to split a supercluster, split the
 * clusters as necessary and return the new cluster indices of this split.
 * Update the end indices of clusters in this supercluster as well.
 */
std::vector<int> split_cluster(
        std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & variant_split_indices,
        std::vector< std::vector<int> > & breakpoints,
        int breakpoint_idx,
        bool print) {
    if (print) printf("Splitting cluster at variant indices (%d, %d)\n", 
            variant_split_indices[0], variant_split_indices[1]);

    std::vector<int> cluster_curr_indices(CALLSETS, 0);
    for (int ci = 0; ci < CALLSETS; ci++) {
        if(print) printf("callset: %s\n", callset_strs[ci].data());

        // find the cluster index corresponding to this variant index
        int var_idx = variant_split_indices[ci];
        auto clust_itr = std::lower_bound(vars[ci]->clusters.begin(),
                vars[ci]->clusters.end(), var_idx);
        int clust_idx = std::distance(vars[ci]->clusters.begin(), clust_itr);
        cluster_curr_indices[ci] = clust_idx;
        if (print) printf("\tvar_idx: %d, clust_idx: %d, *clust_itr: %d\n", var_idx, clust_idx, *clust_itr);

        if (*clust_itr == var_idx) {
            // there is already a cluster break at this variant on this haplotype, do nothing
        } else {
            // we need to split the current cluster in two, save and adjust existing cluster end
            int right_reach = vars[ci]->right_reaches[clust_idx-1];
            int var_pos = vars[ci]->poss[var_idx];
            vars[ci]->right_reaches[clust_idx-1] = var_pos;

            // add new cluster
            vars[ci]->left_reaches.insert(vars[ci]->left_reaches.begin() + clust_idx, var_pos);
            vars[ci]->right_reaches.insert(vars[ci]->right_reaches.begin() + clust_idx, right_reach);
            vars[ci]->clusters.insert(vars[ci]->clusters.begin() + clust_idx, var_idx);
            vars[ci]->nc++;

            // increment breakpoint for all remaining clusters in this supercluster
            for (int bi = breakpoint_idx+1; bi < int(breakpoints.size()); bi++) {
                breakpoints[bi][ci]++;
            }
        }
    }
    return cluster_curr_indices;
}


/******************************************************************************/


/* Get the callset index of the next variant, with a start and end range. */
var_info get_next_variant_info(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & var_curr_indices,
        const std::vector<int> & var_end_indices) {

    int next_callset = -1;
    int next_start_pos = std::numeric_limits<int>::max();
    int next_end_pos = std::numeric_limits<int>::max();
    
    for (int c = 0; c < CALLSETS; c++) {
        if (var_curr_indices[c] < var_end_indices[c]) {
            int start_pos = vars[c]->poss[var_curr_indices[c]];
            int end_pos = vars[c]->poss[var_curr_indices[c]] +
                vars[c]->rlens[var_curr_indices[c]];
            if (start_pos < next_start_pos) {
                next_start_pos = start_pos;
                next_end_pos = end_pos;
                next_callset = c;
            }
        }
    }
    return var_info(next_callset, next_start_pos, next_end_pos);
}


/******************************************************************************/

/* Return the optimal location to split the supercluster.
 */
std::vector<int> get_supercluster_split_location(
        const std::vector< std::shared_ptr<ctgVariants> > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices, bool print) {
    if (print) printf("Finding supercluster split location\n");

    // get original start/end positions of supercluster
    std::vector<int> orig_sc_range = 
        get_supercluster_range(vars, cluster_start_indices, cluster_end_indices);
    int orig_sc_beg_pos = orig_sc_range[0];
    int orig_sc_end_pos = orig_sc_range[1];
    int orig_sc_size = orig_sc_end_pos - orig_sc_beg_pos;

    std::vector<int> var_start_indices(CALLSETS, 0);
    std::vector<int> var_end_indices(CALLSETS, 0);
    for (int c = 0; c < CALLSETS; c++) {
        var_start_indices[c] = vars[c]->clusters[cluster_start_indices[c]];
        var_end_indices[c] = vars[c]->clusters[cluster_end_indices[c]];
    }

    std::vector<int> split_indices = var_start_indices;
    double best_split_score = 0;
    std::vector<int> var_best_split_indices = {};

    // check that there are 2+ variants (this supercluster can be split)
    int total_vars = 0;
    for (int c = 0; c < CALLSETS; c++) {
        total_vars += var_end_indices[c] - var_start_indices[c];
    }
    if (total_vars < 2) return var_best_split_indices; // empty

    // get position and hap of next variant
    var_info curr_var = get_next_variant_info(vars, split_indices, var_end_indices);
    split_indices[curr_var.callset_idx]++;
    var_info next_var = get_next_variant_info(vars, split_indices, var_end_indices);
    while (next_var.callset_idx >= 0) {

        // calculate max split size reduction factor
        int gap = std::max(0, next_var.start_pos - curr_var.end_pos);
        double size_reduction_factor = std::max(double((curr_var.end_pos + gap/2) - orig_sc_beg_pos) / orig_sc_size,
                double(orig_sc_end_pos - (curr_var.end_pos + gap/2)) / orig_sc_size);
        double splits_to_halve_size = -1 / log2(size_reduction_factor);

        // calculate overlap (could weight by number of haps overlapping)
        double split_score = gap / splits_to_halve_size;
        if (print) printf("indices: [%d, %d], gap: %d, frac: %f, splits: %f, score: %f\n",
                split_indices[0], split_indices[1], gap, size_reduction_factor, splits_to_halve_size, split_score);
        if (split_score > best_split_score) {
            best_split_score = split_score;
            var_best_split_indices = split_indices;
        }

        curr_var = next_var;
        split_indices[curr_var.callset_idx]++;
        next_var = get_next_variant_info(vars, split_indices, var_end_indices);
    }
    return var_best_split_indices;
}


/******************************************************************************/


/* Cluster variants using one of two simple heuristic-based methods:
 * 1. "gap N" based clustering
 *     split clusters when genomic distance exceeds N
 *
 * 2. "size N" based clustering
 *     split clusters when genomic distance exceeds max(N, sizeof(X))
 *     for all nearby variants X
 *
 * This function adds `clusters`, `left_reaches`, `right_reaches` metadata
 * to all `VariantData`
 * */
void simple_cluster(std::shared_ptr<variantData> vcf, int callset) {
    bool print = false;

    // cluster each contig
    for (const std::string & ctg : vcf->contigs) {

        // cluster per-haplotype variants: vcf->variants[hap]
        for (int hap = 0; hap < HAPS; hap++) {

            std::shared_ptr<ctgVariants> vars = vcf->variants[hap][ctg];
            if (!vars->n) continue;

            // mark variant boundary between clusters (hence n+1), 
            // with initially one variant per cluster
            std::vector<int> prev_clusters(vars->n+1);
            for (int i = 0; i < vars->n+1; i++) 
                prev_clusters[i] = i;
            std::vector<int> right_reach(vars->n+1), left_reach(vars->n+1);
            std::vector<int> next_clusters, tmp_clusters;

            // save temp clustering (generate_str assumes clustered)
            vars->clusters = prev_clusters;

            // set sentinel reaches
            left_reach[prev_clusters.size()-1] = std::numeric_limits<int>::max();
            right_reach[prev_clusters.size()-1] = std::numeric_limits<int>::max();

            // update all cluster reaches
            for (int var = 0; var < vars->n; var++) {
                int var_size = 0;
                // if g.cluster_method == "gap", keep "var_size" zero
                if (g.cluster_method == "size") {
                    switch (vars->types[var]) {
                        case TYPE_SUB:
                            var_size = 1;
                            break;
                        case TYPE_INS:
                            var_size = vars->alts[var].size();
                            break;
                        case TYPE_DEL:
                            var_size = vars->refs[var].size();
                            break;
                        default:
                            ERROR("Variant type '%s' unexpected at %s:%d.", 
                                    type_strs[vars->types[var]].data(),
                                    ctg.data(), vars->poss[var]);
                    }
                }
                left_reach[var] = vars->poss[var] - 
                        std::max(g.cluster_min_gap, var_size);
                right_reach[var] = vars->poss[var] + vars->rlens[var] +
                        std::max(g.cluster_min_gap, var_size);
                if (print) printf("%s %s:%d hap %d, (%.10s -> %.10s) [%d,%d]\n",
                        callset_strs[callset].data(), ctg.data(), vars->poss[var], hap, 
                        vars->refs[var].data(), vars->alts[var].data(), 
                        left_reach[var], right_reach[var]);
            }

            // merge dependent clusters rightwards
            std::vector<int> tmp_left_reach, tmp_right_reach;
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
                clust += clust_size;
            }
            left_reach.clear();
            right_reach.clear();

            // merge dependent clusters leftwards
            clust = tmp_clusters.size()-1;
            while (clust >= 0) {
                int min_left_reach = tmp_left_reach[clust];
                int max_right_reach = tmp_right_reach[clust];
                while (clust > 0 &&
                        min_left_reach <= tmp_right_reach[clust-1] + g.reach_min_gap) {
                    min_left_reach = std::min(min_left_reach,
                            tmp_left_reach[clust-1]);
                    max_right_reach = std::max(max_right_reach,
                            tmp_right_reach[clust-1]);
                    clust--;
                }
                left_reach.push_back(min_left_reach);
                right_reach.push_back(max_right_reach);
                next_clusters.push_back(tmp_clusters[clust]);
                clust--;
            }
            std::reverse(next_clusters.begin(), next_clusters.end());
            std::reverse(left_reach.begin(), left_reach.end());
            std::reverse(right_reach.begin(), right_reach.end());

            // save final clustering
            vars->clusters = next_clusters;
            vars->nc = int(vars->clusters.size())-1;
            vars->left_reaches = left_reach;
            vars->right_reaches = right_reach;

            if (vars->clusters[vars->nc] != vars->n) {
                ERROR("Mismatch between original and clustered variant count, contig '%s' hap %d",
                        ctg.data(), hap);
            }
            /* printf("CLUSTERS: HAP %d\n", hap); */
            /* for (int i = 0; i < vars->nc; i++) { */
            /*     printf("vars %d-%d: %d-%d\n", vars->clusters[i], vars->clusters[i+1], */
            /*             vars->left_reaches[i], vars->right_reaches[i]); */
            /* } */
        } // hap
    } // ctg
}


/******************************************************************************/


/* Add single-VCF cluster indices to `variantData`. This version assumes that
 * all variant calls are true positives (doesn't allow skipping)
 */
void wf_swg_cluster(variantData * vcf, int ctg_idx, 
        int hap, int sub, int open, int extend) {
    bool print = false;
    std::string ctg = vcf->contigs[ctg_idx];

    // allocate this memory once, use on each cluster
    std::vector<int> offs_buffer(MATS * g.max_size*2 * 
            std::max(sub, open+extend), -2);

    std::shared_ptr<ctgVariants> vars = vcf->variants[hap][ctg];
    if (!vars->n) return;

    // mark variant boundary between clusters (hence n+1), 
    // and whether it was active on the previous iteration.
    // initially, one variant per cluster
    std::vector<int> prev_clusters(vars->n+1);
    for (int i = 0; i < vars->n+1; i++) 
        prev_clusters[i] = i;
    std::vector<bool> prev_active(vars->n+1, true);

    std::vector<int> right_reach(vars->n+1), left_reach(vars->n+1);
    std::vector<int> next_clusters, tmp_clusters;
    std::vector<bool> next_active, tmp_active;

    // while clusters are being merged, loop
    int iter = 0;
    while (std::find(prev_active.begin(), 
                prev_active.end(), true) != prev_active.end()) {
        iter++;
        if (iter > g.max_cluster_itrs) break;
        /* if (print) printf("Iteration %d\n", iter); */

        // count clusters currently being expanded
        int active = 0;
        for (size_t i = 0; i < prev_active.size()-1; i++) {
            if (prev_active[i]) active++;
        }

        // save temp clustering (generate_str assumes clustered)
        vars->clusters = prev_clusters;

        // set sentinel reaches
        left_reach[prev_clusters.size()-1] = std::numeric_limits<int>::max();
        right_reach[prev_clusters.size()-1] = std::numeric_limits<int>::max();

        // update all cluster reaches
        for (size_t clust = 0; clust < prev_clusters.size(); clust++) {

            // only compute necessary reaches (adjacent merge)
            bool left_compute = prev_active[clust];
            bool right_compute = prev_active[clust];

            // no left cluster, or sentinel
            if (clust == prev_clusters.size()-1) {
                left_compute = false;
                right_compute = false;
            }

            // debug print
            if (print && prev_active[clust] && clust < prev_clusters.size()-1) {
                printf("\n\ncluster %d: vars %d-%d, pos %s:%d-%d\n",
                        int(clust), vars->clusters[clust],
                        vars->clusters[clust+1], ctg.data(),
                        vars->poss[vars->clusters[clust]],
                        vars->poss[vars->clusters[clust+1]-1] +
                        vars->rlens[vars->clusters[clust+1]-1]);
                for (int var_idx = vars->clusters[clust]; 
                        var_idx < vars->clusters[clust+1]; var_idx++) {
                    printf("    %s:%d %s %s var:%d\n",
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

            // calculate alignment score
            int l_reach = 0, r_reach = 0, score = 0;
            if (left_compute || right_compute) {

                int beg_idx = vars->clusters[clust];
                int end_idx = vars->clusters[clust+1];
                int beg = std::max(0, vars->poss[beg_idx] - 1);
                int end = std::min(vcf->lengths[ctg_idx], 
                        vars->poss[end_idx-1] + vars->rlens[end_idx-1] + 1);
                std::string query = generate_str(vcf->ref, vars, ctg, beg_idx, end_idx, beg, end);
                std::string ref = vcf->ref->fasta.at(ctg).substr(beg, end-beg);
                wf_swg_align(query, ref, score, sub, open, extend, false);
                if (print) printf("   align score: %d\n", score);
            }

            // LEFT REACH 
            if (left_compute) { // calculate left reach
                std::string query, ref;

                // just after last variant in this cluster
                int beg_pos = vars->poss[vars->clusters[clust]]-1;
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
                int ref_len = end_pos - beg_pos;
                int reach = ref_len - 1;
                while (reach == ref_len-1) { // iterative doubling
                    ref_len *= 2;
                    beg_pos = std::max(0, end_pos - ref_len - std::abs(main_diag) - score/extend - 3); // ensure query is longer
                    query = generate_str(vcf->ref, vars, ctg, 
                                vars->clusters[clust], 
                                vars->clusters[clust+1], 
                                beg_pos, end_pos);
                    ref = vcf->ref->fasta.at(ctg).substr(std::max(0, end_pos-ref_len), ref_len);
                    std::reverse(query.begin(), query.end());
                    std::reverse(ref.begin(), ref.end());
                    // manage buffer for storing offsets
                    size_t offs_size = MATS * (std::max(open+extend, sub)+1) * 
                        (query.size() + ref.size() - 1);
                    if (offs_size > offs_buffer.size())
                        offs_buffer.resize(offs_size, -2);
                    // calculate reach
                    reach = wf_swg_max_reach(query, ref, offs_buffer,
                            main_diag, main_diag_start, score,
                            sub, open, extend, 
                            false /* print */, true /* reverse */);
                    // reset buffer
                    for (size_t i = 0; i < offs_size; i++)
                        offs_buffer[i] = -2;
                    // exit if we've reached the start of the contig
                    if (beg_pos == 0)
                        break;
                }
                if (print) printf("    left reach: %d\n", reach);
                l_reach = end_pos - reach;
                left_reach[clust] = l_reach;

                if (print) {
                    printf("REF:        %s\n", ref.data());
                    printf("QUERY:      %s\n", query.data());
                    printf(" main diag:  %d\n", main_diag);
                    printf("diag start:  %d\n\n", main_diag_start);
                }
            } else {
                l_reach = left_reach[clust];
            }

            // RIGHT REACH
            if (right_compute) { // calculate right reach
                std::string query, ref;

                // right before current cluster
                int beg_pos = vars->poss[vars->clusters[clust]]-1;
                int end_pos = vars->poss[vars->clusters[clust+1]-1] +
                        vars->rlens[vars->clusters[clust+1]-1]+1;

                // get reference end pos of last variant in this cluster
                // (disallow reference main diag matches after this)
                int main_diag_start = vars->poss[ 
                            vars->clusters[clust+1]-1]
                            + vars->rlens[vars->clusters[clust+1]-1] - beg_pos;
                int main_diag = 0;
                for (int vi = vars->clusters[clust];
                        vi < vars->clusters[clust+1]; vi++)
                    main_diag += vars->refs[vi].size() - vars->alts[vi].size();

                // calculate max reaching path to right
                int ref_len = end_pos - beg_pos;
                int reach = ref_len - 1;
                while (reach == ref_len-1) { // iterative doubling
                    ref_len *= 2;
                    end_pos = std::min(vcf->lengths[ctg_idx], 
                            beg_pos + ref_len + std::abs(main_diag) + score/extend + 3);
                    query = generate_str(vcf->ref, vars, ctg, 
                                vars->clusters[clust], 
                                vars->clusters[clust+1], 
                                beg_pos, end_pos);
                    ref = vcf->ref->fasta.at(ctg).substr(beg_pos, std::min(ref_len, end_pos - beg_pos));
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
                    // exit if we've reached the end of the contig
                    if (end_pos == vcf->lengths[ctg_idx])
                        break;
                }
                if (print) printf("   right reach: %d\n", reach);
                r_reach = beg_pos + reach + 1;
                right_reach[clust] = r_reach;

                if (print) {
                    printf("REF:        %s\n", ref.data());
                    printf("QUERY:      %s\n", query.data());
                    printf(" main diag:  %d\n", main_diag);
                    printf("diag start:  %d\n\n", main_diag_start);
                }
            } else {
                r_reach = right_reach[clust];
            }
            if (print) printf("span: %s - %s\n", 
                        l_reach == std::numeric_limits<int>::max() ? 
                            "X" : std::to_string(l_reach).data(), 
                        r_reach == std::numeric_limits<int>::max() ? "X" : std::to_string(r_reach).data());
        }

        // merge dependent clusters rightwards
        std::vector<int> tmp_left_reach, tmp_right_reach;
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
        left_reach.clear();
        right_reach.clear();

        /* printf("tmp clusters:"); */
        /* for (int i = 0; i < 1000; i++) { */
        /*     printf("\t%d", tmp_clusters[i]); */
        /* } */
        /* printf("\n"); */
        /* printf("  tmp active:"); */
        /* for (int i = 0; i < 1000; i++) { */
        /*     printf("\t%d", tmp_active[i] ? 1 : 0); */
        /* } */
        /* printf("\n"); */

        // merge dependent clusters leftwards
        clust = tmp_clusters.size()-1;
        while (clust >= 0) {
            int min_left_reach = tmp_left_reach[clust];
            int max_right_reach = tmp_right_reach[clust];
            bool active = tmp_active[clust];
            while (clust > 0 &&
                    min_left_reach <= tmp_right_reach[clust-1] + g.reach_min_gap) {
                min_left_reach = std::min(min_left_reach,
                        tmp_left_reach[clust-1]);
                max_right_reach = std::max(max_right_reach,
                        tmp_right_reach[clust-1]);
                active = true;
                clust--;
            }
            left_reach.push_back(min_left_reach);
            right_reach.push_back(max_right_reach);
            next_clusters.push_back(tmp_clusters[clust]);
            next_active.push_back(active);
            clust--;
        }
        std::reverse(next_clusters.begin(), next_clusters.end());
        std::reverse(next_active.begin(), next_active.end());
        std::reverse(left_reach.begin(), left_reach.end());
        std::reverse(right_reach.begin(), right_reach.end());

        /* printf("next clusters:"); */
        /* for (int i = 0; i < 1000; i++) { */
        /*     printf("\t%d", next_clusters[i]); */
        /* } */
        /* printf("\n"); */
        /* printf("  next active:"); */
        /* for (int i = 0; i < 1000; i++) { */
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
    }

    // save final clustering
    vars->clusters = prev_clusters;
    vars->nc = int(vars->clusters.size())-1;
    vars->left_reaches = left_reach;
    vars->right_reaches = right_reach;
}
