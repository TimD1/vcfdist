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


/* Merge evaluation results from multiple distinct evaluations. */
void superclusterData::merge(std::shared_ptr<superclusterData> other_data) {
    for (const std::string & ctg : other_data->contigs) {

        // if another contig is completely new, just store the whole thing
        if (std::find(this->contigs.begin(), this->contigs.end(), ctg) == this->contigs.end()) {
            this->superclusters[ctg] = other_data->superclusters[ctg];
            continue;
        }

        // otherwise, contig exists for both, add variants in correct order
        for (int ci = 0; ci < CALLSETS; ci++) {
            int vi1 = 0;
            int vi2 = 0;
            std::shared_ptr<ctgVariants> merged_vars(new ctgVariants(ctg));
            std::shared_ptr<ctgVariants> v1 = this->superclusters[ctg]->callset_vars[ci];
            std::shared_ptr<ctgVariants> v2 = other_data->superclusters[ctg]->callset_vars[ci];
            while (vi1 < v1->n && vi2 < v2->n) {
                if (v1->poss[vi1] < v2->poss[vi2]) { // v1 occurs first
                    merged_vars->add_var(v1, vi1++);
                } else if (v1->poss[vi1] > v2->poss[vi2]) { // v2 occurs first
                    merged_vars->add_var(v2, vi2++);
                } else { // tie, order by variant type
                    // TODO: what if variants are the same ref/alt?
                    if (v1->types[vi1] == TYPE_INS) {
                        merged_vars->add_var(v1, vi1++);
                    } else if (v2->types[vi2] == TYPE_INS) {
                        merged_vars->add_var(v2, vi2++);
                    } else {
                        merged_vars->add_var(v1, vi1++);
                    }
                }
            }
            // add all remaining vars
            while (vi1 < v1->n) { merged_vars->add_var(v1, vi1++); }
            while (vi2 < v2->n) { merged_vars->add_var(v2, vi2++); }
            this->superclusters[ctg]->callset_vars[ci] = merged_vars;
        }
    }
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


// TODO: remove this function after shifting to WFA alignment (with max ED)
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
        int nscs =  std::max(qvars->superclusters[qvars->n-1], tvars->superclusters[tvars->n-1]);

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
            // TODO: this will need to be updated
            // 10 comes from:
            // (4) aln_ptrs    = 1 byte * 2 haps * 2 phasings
            // (2) path_ptrs   = 1 byte * 2 haps
            // (4) path_scores = 2 byte * 2 haps
            // 2 is a fudge factor that I'm adding for now (fragmentation)
            size_t mem = max_lens[QUERY] * max_lens[TRUTH]* 10 * 2;
            double mem_gb = mem / (1000.0 * 1000.0 * 1000.0);
            if (mem_gb > g.max_ram) {
                WARN("Max (%.3fGB) RAM exceeded (%.3fGB req) for supercluster %d, running anyways", 
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
        std::shared_ptr<ctgVariants> > > vars) {

    for (int ctg_idx = 0; ctg_idx < int(this->contigs.size()); ctg_idx++) {
        std::string ctg = contigs[ctg_idx];

        // skip empty contigs
        int nvars = 0;
        for (int h = 0; h < HAPS; h++) nvars += vars[h][ctg]->n;
        if (!nvars) continue;

        // merge variants from each haplotype for this ctg/callset
        std::shared_ptr<ctgVariants> merged_vars(new ctgVariants(ctg));

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
        }

        // initialize indices
        std::vector<int> var_idx(HAPS, 0);
        std::vector<int> clust_idx(HAPS, 0);
        int curr_var_idx = 0;
        int next_left_reach = std::numeric_limits<int>::max();
        int next_right_reach = std::numeric_limits<int>::min();
        while (var_idx[HAP1] < vars[HAP1][ctg]->n ||
                var_idx[HAP2] < vars[HAP2][ctg]->n) {

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
                    /* printf("adding 1|1 var= %s:%d\t%s\t%s\t%s\n", */
                    /*         ctg.data(), */ 
                    /*         vars[HAP1][ctg]->poss[var_idx[HAP1]], */
                    /*         vars[HAP1][ctg]->refs[var_idx[HAP1]].data(), */
                    /*         vars[HAP1][ctg]->alts[var_idx[HAP1]].data(), */
                    /*         gt_strs[vars[HAP1][ctg]->orig_gts[var_idx[HAP1]]].data() */
                    /*         ); */
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
                    /* printf("adding %s var= %s:%d\t%s\t%s\t%s\n", */
                    /*         hap_idx == 0 ? "1|0" : "0|1", */
                    /*         ctg.data(), */ 
                    /*         vars[hap_idx][ctg]->poss[var_idx[hap_idx]], */
                    /*         vars[hap_idx][ctg]->refs[var_idx[hap_idx]].data(), */
                    /*         vars[hap_idx][ctg]->alts[var_idx[hap_idx]].data(), */
                    /*         gt_strs[vars[hap_idx][ctg]->orig_gts[var_idx[hap_idx]]].data() */
                    /*             ); */
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
                    /* printf("adding 1|0 var= %s:%d\t%s\t%s\t%s\n", */
                    /*         ctg.data(), */ 
                    /*         vars[HAP1][ctg]->poss[var_idx[HAP1]], */
                    /*         vars[HAP1][ctg]->refs[var_idx[HAP1]].data(), */
                    /*         vars[HAP1][ctg]->alts[var_idx[HAP1]].data(), */
                    /*         gt_strs[vars[HAP1][ctg]->orig_gts[var_idx[HAP1]]].data() */
                    /*             ); */
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
                    /* printf("adding 0|1 var= %s:%d\t%s\t%s\t%s\n", */
                    /*         ctg.data(), */ 
                    /*         vars[HAP2][ctg]->poss[var_idx[HAP2]], */
                    /*         vars[HAP2][ctg]->refs[var_idx[HAP2]].data(), */
                    /*         vars[HAP2][ctg]->alts[var_idx[HAP2]].data(), */
                    /*         gt_strs[vars[HAP1][ctg]->orig_gts[var_idx[HAP1]]].data() */
                    /*             ); */
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
                    if (vars[h][ctg]->left_reaches[clust_idx[h]] <= 
                            vars[h^1][ctg]->left_reaches[clust_idx[h^1]]) {
                        next_left_reach = vars[h][ctg]->left_reaches[clust_idx[h]];
                        next_right_reach = std::max(next_right_reach,
                            vars[h][ctg]->right_reaches[clust_idx[h]]);
                    }
                }
            }
            /* printf("curr = (%d, %d)\tnext = (%d, %d)\n", */
            /*         curr_left_reach, curr_right_reach, */
            /*         next_left_reach, next_right_reach); */

            // starting new supercluster
            if (next_left_reach > curr_right_reach) {
                /* printf("new supercluster, adding curr = %d (%d, %d)\n\n", */
                /*         curr_var_idx, curr_left_reach, curr_right_reach); */

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
                            vars[h][ctg]->left_reaches[clust_idx[h]] <=
                            vars[h^1][ctg]->left_reaches[clust_idx[h^1]]) {
                        next_left_reach  = std::min(next_left_reach,
                                vars[h][ctg]->left_reaches[clust_idx[h]]);
                        next_right_reach = std::max(next_right_reach,
                                vars[h][ctg]->right_reaches[clust_idx[h]]);
                    }
                }

            } else { // same supercluster
                /* printf("same supercluster\n"); */
                curr_left_reach = std::min(curr_left_reach, next_left_reach);
                curr_right_reach = std::max(curr_right_reach, next_right_reach);
            }

        } // while variants remain

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

    // supercluster
    this->supercluster();
}


/******************************************************************************/

/* Supercluster using left_reach and right_reach of each cluster, calculated
 * during wf_swg_cluster().
 */
void superclusterData::supercluster() {

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
        auto const & vars = this->superclusters[ctg]->callset_vars;
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

            // get supercluster start/end positions (allowing empty haps)
            int beg_pos = std::numeric_limits<int>::max();
            int end_pos = -1;
            for (int c = 0; c < CALLSETS; c++) {
                if (next_brks[c] - brks[c]) {
                    beg_pos = std::min(beg_pos,
                        vars[c]->poss[
                            vars[c]->clusters[brks[c]]]-1);
                    end_pos = std::max(end_pos,
                            vars[c]->poss[
                                vars[c]->clusters[next_brks[c]]-1] +
                            vars[c]->rlens[
                                vars[c]->clusters[next_brks[c]]-1] + 1);
                }
            }

            // update summary metrics
            largest_supercluster = std::max(largest_supercluster, end_pos-beg_pos);
            total_bases += end_pos-beg_pos;
            int this_vars = 0;
            for (int c = 0; c < CALLSETS; c++) {
                if (vars[c]->nc)
                    this_vars += vars[c]->clusters[next_brks[c]] - vars[c]->clusters[brks[c]];
            }
            most_vars = std::max(most_vars, this_vars);
            total_vars += this_vars;

            // debug print
            if (false) {
                printf("\nSUPERCLUSTER: %d\n", sc_idx);
                printf("POS: %s:%d-%d\n", ctg.data(), beg_pos, end_pos);
                printf("SIZE: %d\n", end_pos - beg_pos);
                for (int c = 0; c < CALLSETS; c++) {
                    printf("%s: clusters %d-%d of %d\n",
                            callset_strs[c].data(),
                            brks[c], next_brks[c], vars[c]->nc);
                    if (vars[c]->nc)
                    for (int v = vars[c]->clusters[ brks[c]]; 
                            v < vars[c]->clusters[ next_brks[c] ]; v++) {
                        printf("    var %d = %s:%d\t%s\t%s\t%s\n", v, ctg.data(),
                                vars[c]->poss[v], 
                                vars[c]->refs[v].data(),
                                vars[c]->alts[v].data(),
                                gt_strs[vars[c]->orig_gts[v]].data());
                    }
                }
                printf("curr_right: %d, lefts: %d %d\n", curr_right, lefts[0], lefts[1]);
                printf("variants: %d, bases: %d\n", this_vars, end_pos-beg_pos);
            }

            // save supercluster information
            for (int ci = 0; ci < CALLSETS; ci++) {
                for (int vi = vars[ci]->clusters[brks[ci]]; 
                         vi < vars[ci]->clusters[next_brks[ci]]; vi++) {
                    vars[ci]->superclusters[vi] = sc_idx;
                }
            }
            sc_idx++;

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
    for (std::string ctg : vcf->contigs) {

        // cluster per-haplotype variants: vcf->variants[hap]
        for (int hap = 0; hap < HAPS; hap++) {

            std::shared_ptr<ctgVariants> vars = vcf->variants[hap][ctg];
            if (!vars->n) break;

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
                std::vector< std::vector< std::vector<uint8_t> > > ptrs(MATS);
                std::vector< std::vector< std::vector<int> > > offs(MATS);
                wf_swg_align(query, ref, ptrs, offs, score, sub, open, extend, false);
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
                    ref = vcf->ref->fasta.at(ctg).substr(end_pos-ref_len, ref_len);
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
                    ref = vcf->ref->fasta.at(ctg).substr(beg_pos, ref_len);
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
