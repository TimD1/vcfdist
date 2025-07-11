#include <algorithm>
#include <cmath>
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
    this->sc_phase.push_back(PHASE_NONE);
    this->pb_phase.push_back(PHASE_NONE);
    this->orig_phase_dist.push_back(-1);
    this->swap_phase_dist.push_back(-1);
    this->n++;
}

void ctgSuperclusters::set_phase(
           int sc_idx, int phase, 
           int orig_phase_dist, 
           int swap_phase_dist) {
    this->sc_phase[sc_idx] = phase;
    this->orig_phase_dist[sc_idx] = orig_phase_dist;
    this->swap_phase_dist[sc_idx] = swap_phase_dist;
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
        for (int sc_idx = 0; sc_idx < ctg_scs->n; sc_idx++) {
            int max_query_len = 0;
            // get max query length
            for (int hap = 0; hap < HAPS; hap++) {
                // get range of included variants
                /* printf("QUERY hap=%d\n", hap); */
                if (ctg_scs->ctg_variants[QUERY][hap]->clusters.size() == 0) continue;
                /* printf("nscs=%d, sc_size=%d, sc_idx=%d\n", ctg_scs->n, int(scs.size()), sc_idx); */
                /* printf("clust_beg=%d, clust_end=%d, clust_size=%d\n", scs[sc_idx], scs[sc_idx+1], int(vars->clusters.size())); */
                int var_beg = ctg_scs->ctg_variants[QUERY][hap]->clusters[
                        ctg_scs->superclusters[QUERY][hap][sc_idx]];
                int var_end = ctg_scs->ctg_variants[QUERY][hap]->clusters[
                        ctg_scs->superclusters[QUERY][hap][sc_idx+1]];
                /* printf("var_beg=%d, var_end=%d, nvar=%d, var_size=%d\n", var_beg, var_end, vars->n, int(vars->poss.size())); */
                // calculate query len after applying variants
                int query_len = ctg_scs->ends[sc_idx] - ctg_scs->begs[sc_idx];
                for (int var = var_beg; var < var_end; var++) {
                    query_len += (ctg_scs->ctg_variants[QUERY][hap]->alts[var].size() - 
                            ctg_scs->ctg_variants[QUERY][hap]->refs[var].size());
                }
                max_query_len = std::max(max_query_len, query_len);
            }
            
            // get max truth length
            int max_truth_len = 0;
            for (int hap = 0; hap < HAPS; hap++) {
                // get range of included variants
                if (ctg_scs->ctg_variants[TRUTH][hap]->clusters.size() == 0) continue;
                int var_beg = ctg_scs->ctg_variants[TRUTH][hap]->clusters[
                        ctg_scs->superclusters[TRUTH][hap][sc_idx]];
                int var_end = ctg_scs->ctg_variants[TRUTH][hap]->clusters[
                        ctg_scs->superclusters[TRUTH][hap][sc_idx+1]];
                // calculate truth len after applying variants
                int truth_len = ctg_scs->ends[sc_idx] - ctg_scs->begs[sc_idx];
                for (int var = var_beg; var < var_end; var++) {
                    truth_len += (ctg_scs->ctg_variants[TRUTH][hap]->alts[var].size() - 
                            ctg_scs->ctg_variants[TRUTH][hap]->refs[var].size());
                }
                max_truth_len = std::max(max_truth_len, truth_len);
            }

            // calculate memory usage
            // 10 comes from:
            // (4) aln_ptrs    = 1 byte * 2 haps * 2 phasings
            // (2) path_ptrs   = 1 byte * 2 haps
            // (4) path_scores = 2 byte * 2 haps
            // 2 is a fudge factor that I'm adding for now (fragmentation)
            size_t mem = size_t(max_query_len) * size_t(max_truth_len) * 10 * 2;
            double mem_gb = mem / (1000.0 * 1000.0 * 1000.0);
            if (mem_gb > g.max_ram) {
                WARN("Max (%.3fGB) RAM exceeded (%.3fGB req) for supercluster %s:%d-%d, running anyways. Lower --max-supercluster-size if this fails.", 
                        g.max_ram, mem_gb, ctg.data(), ctg_scs->begs[sc_idx], ctg_scs->ends[sc_idx]);
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

superclusterData::superclusterData(
        std::shared_ptr<variantData> query_ptr,
        std::shared_ptr<variantData> truth_ptr,
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

    // set pointers to variant lists (per contig)
    for (std::string ctg : this->contigs) {
        try {
            this->superclusters[ctg]->ctg_variants[QUERY][HAP1] = 
                    query_ptr->variants[HAP1][ctg];
            this->superclusters[ctg]->ctg_variants[QUERY][HAP2] = 
                    query_ptr->variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Query VCF does not contain contig '%s'", ctg.data());
        }
        try {
            this->superclusters[ctg]->ctg_variants[TRUTH][HAP1] = 
                    truth_ptr->variants[HAP1][ctg];
            this->superclusters[ctg]->ctg_variants[TRUTH][HAP2] = 
                    truth_ptr->variants[HAP2][ctg];
        } catch (const std::exception & e) {
            ERROR("Truth VCF does not contain contig '%s'", ctg.data());
        }
    }
    this->supercluster(false);
    this->transfer_phase_sets();
}


/******************************************************************************/

/* Using per-variant phase set information from PS tags, calculate per-supercluster
 * phase sets. */
void superclusterData::transfer_phase_sets() {
    bool print = false;

    int total_query_phase_sets = 0;
    int total_truth_phase_sets = 0;
    int total_multiple_phase_sets = 0;
    int total_non_increasing = 0;
    std::vector<int> query_phase_set_sizes;
    std::vector<int> truth_phase_set_sizes;
    int ctg_idx = 0;
    for (std::string ctg : this->contigs) { // for each contig
        std::shared_ptr<ctgSuperclusters> ctg_scs = this->superclusters[ctg];
        if (print) printf("ctg: %s\n", ctg.data());
        if (!ctg_scs->n) continue;

        // set convenience variables
        int q1i = 0; int q2i = 0; int t1i = 0; int t2i = 0; // indices
        std::shared_ptr<ctgVariants> q1v = ctg_scs->ctg_variants[QUERY][HAP1]; // vars
        std::shared_ptr<ctgVariants> q2v = ctg_scs->ctg_variants[QUERY][HAP2];
        std::shared_ptr<ctgVariants> t1v = ctg_scs->ctg_variants[TRUTH][HAP1];
        std::shared_ptr<ctgVariants> t2v = ctg_scs->ctg_variants[TRUTH][HAP2];
        std::vector<int> & q1sc = ctg_scs->superclusters[QUERY][HAP1]; // superclusters
        std::vector<int> & q2sc = ctg_scs->superclusters[QUERY][HAP2];
        std::vector<int> & t1sc = ctg_scs->superclusters[TRUTH][HAP1];
        std::vector<int> & t2sc = ctg_scs->superclusters[TRUTH][HAP2];

        // get first phase set (to backfill all preceding zeros)
        int first_pos = std::numeric_limits<int>::max();
        int phase_set = 0;
        for (; q1i < q1v->n; q1i++) {
            if (q1v->phase_sets[q1i] != 0) {
                if (q1v->poss[q1i] < first_pos) {
                    first_pos = q1v->poss[q1i];
                    phase_set = q1v->phase_sets[q1i];
                }
                break;
            }
        }
        for (; q2i < q2v->n; q2i++) {
            if (q2v->phase_sets[q2i] != 0) {
                if (q2v->poss[q2i] < first_pos) {
                    first_pos = q2v->poss[q2i];
                    phase_set = q2v->phase_sets[q2i];
                }
                break;
            }
        }
        for (; t1i < t1v->n; t1i++) {
            if (t1v->phase_sets[t1i] != 0) {
                if (t1v->poss[t1i] < first_pos) {
                    first_pos = t1v->poss[t1i];
                    phase_set = t1v->phase_sets[t1i];
                }
                break;
            }
        }
        for (; t2i < t2v->n; t2i++) {
            if (t2v->phase_sets[t2i] != 0) {
                if (t2v->poss[t2i] < first_pos) {
                    first_pos = t2v->poss[t2i];
                    phase_set = t2v->phase_sets[t2i];
                }
                break;
            }
        }
        if (print) printf("first phase set: %d\n", phase_set);

        // carry phase_set across superclusters (if no heterozygous variants)
        int query_phase_set = 0;
        int truth_phase_set = 0;
        int query_ps_beg = 0; int query_ps_end = 0;
        int truth_ps_beg = 0; int truth_ps_end = 0;
        for (int sci = 0; sci < ctg_scs->n; sci++) { // for each supercluster
            int query_ps_ct = 0;
            int truth_ps_ct = 0;
            int non_increasing = 0;
            if (print) printf("supercluster: %d\n", sci);

            // check all variants in supercluster for each haplotype
            if (q1v->n)
            for (q1i = q1v->clusters[q1sc[sci]]; // QUERY HAP 1
                    q1i < q1v->clusters[q1sc[sci+1]]; q1i++) {
                if (q1v->phase_sets[q1i]) { // non-zero, has PS tag
                    if (q1v->phase_sets[q1i] > query_phase_set) { // new PS
                        if (query_phase_set) query_phase_set_sizes.push_back(query_ps_end - query_ps_beg);
                        phase_set = query_phase_set = q1v->phase_sets[q1i];
                        total_query_phase_sets++; query_ps_ct++;
                        if (print) printf(" Q1 PS:%d %d-%d", phase_set, query_ps_beg, query_ps_end);
                        query_ps_beg = q1v->poss[q1i];
                        query_ps_end = q1v->poss[q1i] + q1v->rlens[q1i];
                    } else if (q1v->phase_sets[q1i] == query_phase_set) { // same 
                        if (!query_ps_ct) { query_ps_ct++; }
                        query_ps_end = std::max(query_ps_end, q1v->poss[q1i] + q1v->rlens[q1i]);
                    } else {
                        query_ps_ct++;
                        non_increasing++;
                    }
                }
            }
            if (q2v->n)
            for (q2i = q2v->clusters[q2sc[sci]]; // QUERY HAP 2
                    q2i < q2v->clusters[q2sc[sci+1]]; q2i++) {
                if (q2v->phase_sets[q2i]) { // non-zero, has PS tag
                    if (q2v->phase_sets[q2i] > query_phase_set) { // new
                        if (query_phase_set) query_phase_set_sizes.push_back(query_ps_end - query_ps_beg);
                        phase_set = query_phase_set = q2v->phase_sets[q2i];
                        total_query_phase_sets++; query_ps_ct++;
                        if (print) printf(" Q2 PS:%d %d-%d", phase_set, query_ps_beg, query_ps_end);
                        query_ps_beg = q2v->poss[q2i];
                        query_ps_end = q2v->poss[q2i] + q2v->rlens[q2i];
                    } else if (q2v->phase_sets[q2i] == query_phase_set) { // same 
                        if (!query_ps_ct) { query_ps_ct++; }
                        query_ps_end = std::max(query_ps_end, q2v->poss[q2i] + q2v->rlens[q2i]);
                    } else {
                        query_ps_ct++;
                        non_increasing++;
                    }
                }
            }
            if (t1v->n)
            for (t1i = t1v->clusters[t1sc[sci]]; // TRUTH HAP 1
                    t1i < t1v->clusters[t1sc[sci+1]]; t1i++) {
                if (t1v->phase_sets[t1i]) { // non-zero, has PS tag
                    if (t1v->phase_sets[t1i] > truth_phase_set) { // new
                        if (truth_phase_set) truth_phase_set_sizes.push_back(truth_ps_end - truth_ps_beg);
                        phase_set = truth_phase_set = t1v->phase_sets[t1i];
                        total_truth_phase_sets++; truth_ps_ct++;
                        if (print) printf(" T1 PS:%d %d-%d", phase_set, truth_ps_beg, truth_ps_end);
                        truth_ps_beg = t1v->poss[t1i];
                        truth_ps_end = t1v->poss[t1i] + t1v->rlens[t1i];
                    } else if (t1v->phase_sets[t1i] == truth_phase_set) { // same 
                        if (!truth_ps_ct) { truth_ps_ct++; }
                        truth_ps_end = std::max(truth_ps_end, t1v->poss[t1i] + t1v->rlens[t1i]);
                    } else {
                        truth_ps_ct++;
                        non_increasing++;
                    }
                }
            }
            if (t2v->n)
            for (t2i = t2v->clusters[t2sc[sci]]; // TRUTH HAP 2
                    t2i < t2v->clusters[t2sc[sci+1]]; t2i++) {
                if (t2v->phase_sets[t2i]) { // non-zero, has PS tag
                    if (t2v->phase_sets[t2i] > truth_phase_set) { // new
                        if (truth_phase_set) truth_phase_set_sizes.push_back(truth_ps_end - truth_ps_beg);
                        phase_set = truth_phase_set = t2v->phase_sets[t2i];
                        total_truth_phase_sets++; truth_ps_ct++;
                        if (print) printf(" T2 PS:%d %d-%d", phase_set, truth_ps_beg, truth_ps_end);
                        truth_ps_beg = t2v->poss[t2i];
                        truth_ps_end = t2v->poss[t2i] + t2v->rlens[t2i];
                    } else if (t2v->phase_sets[t2i] == truth_phase_set) { // same 
                        if (!truth_ps_ct) { truth_ps_ct++; }
                        truth_ps_end = std::max(truth_ps_end, t2v->poss[t2i] + t2v->rlens[t2i]);
                    } else {
                        truth_ps_ct++;
                        non_increasing++;
                    }
                }
            }

            // if one supercluster contains variants in multiple phase sets
            if (query_ps_ct > 1 || truth_ps_ct > 1) {
                total_multiple_phase_sets++;
            } else {
                total_non_increasing += non_increasing;
            }

            // add phase_set to supercluster
            ctg_scs->phase_sets.push_back(phase_set);
            if (print) printf(" final=%d\n", phase_set);
        }

        // add final phase set on contig
        if (!query_phase_set) {
            query_phase_set_sizes.push_back(this->lengths[ctg_idx]);
            total_query_phase_sets++;
        } else {
            query_phase_set_sizes.push_back(query_ps_end - query_ps_beg);
            total_query_phase_sets++;
        }
        if (!truth_phase_set) {
            truth_phase_set_sizes.push_back(this->lengths[ctg_idx]);
            total_truth_phase_sets++;
        } else {
            truth_phase_set_sizes.push_back(truth_ps_end - truth_ps_beg);
            total_truth_phase_sets++;
        }
        ctg_idx++;
    }

    // calculate phaseset NG50
    size_t total_bases = 0;
    for (size_t i = 0; i < this->contigs.size(); i++) {
        total_bases += lengths[i];
    }

    int query_pb_ng50 = calc_ng50(query_phase_set_sizes, total_bases);
    int truth_pb_ng50 = calc_ng50(truth_phase_set_sizes, total_bases);

    // print
    if (total_multiple_phase_sets)
        WARN("%d total superclusters contain variants from multiple phase sets (PS)",
                total_multiple_phase_sets)

    if (total_non_increasing)
        WARN("%d total phase sets (PS) are in non-increasing order",
                total_non_increasing)

    if (g.verbosity >= 1) INFO("              QUERY phase sets: %d", total_query_phase_sets);
    if (g.verbosity >= 1) INFO("        QUERY phase block NG50: %d", query_pb_ng50);
    if (g.verbosity >= 1) INFO("              TRUTH phase sets: %d", total_truth_phase_sets);
    if (g.verbosity >= 1) INFO("        TRUTH phase block NG50: %d", truth_pb_ng50);
}


/* Supercluster using left_reach and right_reach of each cluster, calculated
 * during wf_swg_cluster().
 */
void superclusterData::supercluster(bool print) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[4/8] Superclustering TRUTH and QUERY variants%s",
            COLOR_PURPLE, COLOR_WHITE);

    // iterate over each contig
    int total_superclusters = 0;
    int largest_supercluster = 0;
    int total_vars = 0;
    int most_vars = 0;
    int total_bases = 0;
    for (std::string ctg : this->contigs) {

        // skip empty contigs
        int nvars = 0;
        for (int i = 0; i < CALLSETS*HAPS; i++) {
            nvars += this->superclusters[ctg]->ctg_variants[i>>1][i&1]->n;
        }
        if (!nvars) continue;

        // for each cluster of variants (merge query and truth haps)
        auto & vars = this->superclusters[ctg]->ctg_variants;
        std::vector<int> brks = {0, 0, 0, 0}; // start of current supercluster
        while (true) {

            // init: empty supercluster
            std::vector<int> next_brks = brks; // end of current supercluster
            std::vector<int> lefts(4, std::numeric_limits<int>::max());
            for (int i = 0; i < CALLSETS*HAPS; i++) {
                if (brks[i] < int(vars[i>>1][i&1]->clusters.size())-1) {
                    lefts[i] = vars[i>>1][i&1]->left_reaches[ next_brks[i] ];
                }
            }

            // get first cluster, end if all haps are off end
            int idx = std::distance(lefts.begin(),
                    std::min_element(lefts.begin(), lefts.end()));
            if (lefts[idx] == std::numeric_limits<int>::max()) break;

            // initialize cluster merging with first to start
            int curr_right = vars[idx>>1][idx&1]->right_reaches[ next_brks[idx] ];
            next_brks[idx]++;
            lefts[idx] = next_brks[idx] < int(vars[idx>>1][idx&1]->clusters.size())-1 ?
                    vars[idx>>1][idx&1]->left_reaches[ next_brks[idx] ] :
                    std::numeric_limits<int>::max();

            // keep expanding cluster while possible
            bool just_active = true;
            while (just_active) {
                just_active = false;
                for (int i = 0; i < CALLSETS*HAPS; i++) {
                    while (lefts[i] <= curr_right) {
                        curr_right = std::max(curr_right,
                            vars[i>>1][i&1]->right_reaches[ next_brks[i] ]);
                        next_brks[i]++;
                        lefts[i] = next_brks[i] < int(vars[i>>1][i&1]->clusters.size())-1 ?
                            vars[i>>1][i&1]->left_reaches[ next_brks[i] ] :
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
                printf("\nSUPERCLUSTER: %d\n", this->superclusters[ctg]->n);
                printf("POS: %s:%d-%d\n", ctg.data(), beg_pos, end_pos);
                printf("SIZE: %d\n", end_pos - beg_pos);
                for (int i = 0; i < CALLSETS*HAPS; i++) {
                    printf("%s%d: clusters %d-%d of %d\n",
                            callset_strs[i>>1].data(), (i&1)+1,
                            brks[i], next_brks[i], 
                            std::max(int(vars[i>>1][i&1]->clusters.size())-1, 0));
                    if (vars[i>>1][i&1]->clusters.size())
                    for (int v = vars[i>>1][i&1]->clusters[ brks[i]]; v < vars[i>>1][i&1]->clusters[ next_brks[i] ]; v++) {
                        printf("    var %d = %s:%d\t%s\t%s\n", v, ctg.data(),
                                vars[i>>1][i&1]->poss[v], 
                                vars[i>>1][i&1]->refs[v].data(),
                                vars[i>>1][i&1]->alts[v].data());
                    }
                }
                printf("curr_right: %d, lefts: %d %d %d %d\n",
                    curr_right, lefts[0], lefts[1], lefts[2], lefts[3]);
            }

            // split large supercluster if necessary
            if (end_pos - beg_pos > g.max_supercluster_size) {
                std::vector< std::vector<int> > all_brks = split_large_supercluster(vars, brks, next_brks, print);
                WARN("Max supercluster size (%d) exceeded (%d) at %s:%d-%d, breaking up into %d superclusters", 
                        g.max_supercluster_size, end_pos-beg_pos, ctg.data(), beg_pos, end_pos, 
                        int(all_brks.size()-1));
                for (int brk_idx = 0; brk_idx < int(all_brks.size())-1; brk_idx++) {
                    brks = all_brks[brk_idx];
                    next_brks = all_brks[brk_idx+1];
                    poss = get_supercluster_range(vars, brks, next_brks);
                    beg_pos = poss[0];
                    end_pos = poss[1];
                    this->superclusters[ctg]->add_supercluster(brks, beg_pos, end_pos);

                    // update summary metrics
                    largest_supercluster = std::max(largest_supercluster, end_pos-beg_pos);
                    total_bases += end_pos-beg_pos;
                    int this_vars = 0;
                    for (int i = 0; i < CALLSETS*HAPS; i++) {
                        if (vars[i>>1][i&1]->clusters.size())
                            this_vars += vars[i>>1][i&1]->clusters[ next_brks[i] ] - 
                                vars[i>>1][i&1]->clusters[ brks[i] ];
                    }
                    most_vars = std::max(most_vars, this_vars);
                    total_vars += this_vars;
                }

            } else {
                // save alignment information
                this->superclusters[ctg]->add_supercluster(brks, beg_pos, end_pos);

				// update summary statistics
				total_bases += end_pos-beg_pos;
				int this_vars = 0;
				for (int i = 0; i < CALLSETS*HAPS; i++) {
					if (vars[i>>1][i&1]->clusters.size())
						this_vars += vars[i>>1][i&1]->clusters[ next_brks[i] ] - 
							vars[i>>1][i&1]->clusters[ brks[i] ];
				}
				most_vars = std::max(most_vars, this_vars);
				total_vars += this_vars;
            }

            // reset for next active cluster
            brks = next_brks;
        }

        // add sentinel, not actually a supercluster
        this->superclusters[ctg]->add_supercluster(brks, 
            std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
        this->superclusters[ctg]->n--;

        total_superclusters += this->superclusters[ctg]->n;
    }
    if (g.verbosity >= 1) INFO("           Total superclusters: %d", total_superclusters);
    if (g.verbosity >= 1) INFO("  Largest supercluster (bases): %d", largest_supercluster);
    if (g.verbosity >= 1) INFO("  Largest supercluster  (vars): %d", most_vars);
    if (g.verbosity >= 1 && total_superclusters) INFO("  Average supercluster (bases): %.3f", total_bases / float(total_superclusters));
    if (g.verbosity >= 1 && total_superclusters) INFO("  Average supercluster  (vars): %.3f", total_vars / float(total_superclusters));
}


/******************************************************************************/


/* Calculate the size of a supercluster from a list of variants and the start/stop indices of the
 * clusters that compose the supercluster (on each hap).
 * Returns a vector containing the (start, end) positions.
 */
std::vector<int> get_supercluster_range(
        const std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & cluster_start_indices, // inclusive
        const std::vector<int> & cluster_end_indices) { // exclusive

    int beg_pos = std::numeric_limits<int>::max();
    int end_pos = -1;

    for (int i = 0; i < CALLSETS*HAPS; i++) {
        // if there is a cluster on this hap, update beginning and end positions
        if (cluster_end_indices[i] - cluster_start_indices[i]) {
            if (cluster_start_indices[i] >= int(vars[i>>1][i&1]->clusters.size())) {
                ERROR("Cluster start indices invalid in get_supercluster_range()");
            }
            if (cluster_end_indices[i] >= int(vars[i>>1][i&1]->clusters.size())) {
                ERROR("Cluster end indices invalid in get_supercluster_range()");
            }

            // one position left of the leftmost variant
            beg_pos = std::min(beg_pos,
                vars[i>>1][i&1]->poss[
                    vars[i>>1][i&1]->clusters[cluster_start_indices[i]]]-1);

            // one position right of the rightmost variant end in the last included cluster
            end_pos = std::max(end_pos,
                    vars[i>>1][i&1]->poss[
                        vars[i>>1][i&1]->clusters[cluster_end_indices[i]]-1] +
                    vars[i>>1][i&1]->rlens[
                        vars[i>>1][i&1]->clusters[cluster_end_indices[i]]-1] + 1);
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
        std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
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

                if (int(best_var_split.size()) == HAPS*CALLSETS) { // found valid split
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

/* Given a 4-tuple of variant indices (on each hap), at which to split a supercluster, split the
 * clusters as necessary and return the new cluster indices of this split.
 * Update the end indices of clusters in this supercluster as well.
 */
std::vector<int> split_cluster(
        std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & variant_split_indices,
        std::vector< std::vector<int> > & breakpoints,
        int breakpoint_idx,
        bool print) {
    if (print) printf("Splitting cluster at (%d, %d, %d, %d)\n", 
            variant_split_indices[0], variant_split_indices[1],
            variant_split_indices[2], variant_split_indices[3]);

    std::vector<int> cluster_curr_indices(CALLSETS*HAPS, 0);
    for (int i = 0; i < CALLSETS*HAPS; i++) {
        if(print) printf("index: %d\n", i);
        if (not vars[i>>1][i&1]->clusters.size()) continue;

        // find the cluster index corresponding to this variant index
        int var_idx = variant_split_indices[i];
        auto clust_itr = std::lower_bound(vars[i>>1][i&1]->clusters.begin(),
                vars[i>>1][i&1]->clusters.end(), var_idx);
        int clust_idx = std::distance(vars[i>>1][i&1]->clusters.begin(), clust_itr);
        cluster_curr_indices[i] = clust_idx;
        if (print) printf("\tvar_idx: %d, clust_idx: %d, *clust_itr: %d\n", var_idx, clust_idx, *clust_itr);

        if (*clust_itr == var_idx) {
            // there is already a cluster break at this variant on this haplotype, do nothing
        } else {
            // we need to split the current cluster in two, save and adjust existing cluster end
            int right_reach = vars[i>>1][i&1]->right_reaches[clust_idx-1];
            int var_pos = vars[i>>1][i&1]->poss[var_idx];
            vars[i>>1][i&1]->right_reaches[clust_idx-1] = var_pos;

            // add new cluster
            vars[i>>1][i&1]->left_reaches.insert(vars[i>>1][i&1]->left_reaches.begin() + clust_idx, var_pos);
            vars[i>>1][i&1]->right_reaches.insert(vars[i>>1][i&1]->right_reaches.begin() + clust_idx,
                    right_reach);
            vars[i>>1][i&1]->clusters.insert(vars[i>>1][i&1]->clusters.begin() + clust_idx, var_idx);

            // increment breakpoint for all remaining clusters in this supercluster
            for (int j = breakpoint_idx+1; j < int(breakpoints.size()); j++) {
                breakpoints[j][i]++;
            }
        }
    }
    return cluster_curr_indices;
}


/******************************************************************************/


/* Get the index of the next variant, with a start and end range defined for each haplotype.
   NOTE: this is currently always 4: truth and query, hap1 and hap2
 */
var_info get_next_variant_info(
        const std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & var_curr_indices,
        const std::vector<int> & var_end_indices) {

    int next_var_idx = -1;
    int next_start_pos = std::numeric_limits<int>::max();
    int next_end_pos = std::numeric_limits<int>::max();
    
    for (int i = 0; i < CALLSETS*HAPS; i++) {
        if (var_curr_indices[i] < var_end_indices[i]) {
            int start_pos = vars[i>>1][i&1]->poss[var_curr_indices[i]];
            int end_pos = vars[i>>1][i&1]->poss[var_curr_indices[i]] +
                vars[i>>1][i&1]->rlens[var_curr_indices[i]];
            if (start_pos < next_start_pos) {
                next_start_pos = start_pos;
                next_end_pos = end_pos;
                next_var_idx = i;
            }
        }
    }
    return var_info(next_var_idx, next_start_pos, next_end_pos);
}


/******************************************************************************/

/* Return the optimal location to split the supercluster.
 */
std::vector<int> get_supercluster_split_location(
        const std::vector< std::vector< std::shared_ptr<ctgVariants> > > & vars,
        const std::vector<int> & cluster_start_indices,
        const std::vector<int> & cluster_end_indices, bool print) {
    if (print) printf("Finding supercluster split location\n");

    // get original start/end positions of supercluster
    std::vector<int> orig_sc_range = 
        get_supercluster_range(vars, cluster_start_indices, cluster_end_indices);
    int orig_sc_beg_pos = orig_sc_range[0];
    int orig_sc_end_pos = orig_sc_range[1];
    int orig_sc_size = orig_sc_end_pos - orig_sc_beg_pos;
    if (print) printf("original range: %d-%d, size %d\n", orig_sc_beg_pos, orig_sc_end_pos, orig_sc_size);

    std::vector<int> var_start_indices(CALLSETS*HAPS, 0);
    std::vector<int> var_end_indices(CALLSETS*HAPS, 0);
    for (int i = 0; i < CALLSETS*HAPS; i++) {
        if (vars[i>>1][i&1]->clusters.size()) {
            var_start_indices[i] = vars[i>>1][i&1]->clusters[cluster_start_indices[i]];
            var_end_indices[i] = vars[i>>1][i&1]->clusters[cluster_end_indices[i]];
            if (print) printf("start/end var indices on %s hap %d: %d-%d\n",
                    callset_strs[i>>1].data(), i&1, var_start_indices[i], var_end_indices[i]);
        }
    }

    std::vector<int> split_indices = var_start_indices;
    double best_split_score = 0;
    std::vector<int> var_best_split_indices = {};

    // check that there are 2+ variants (this supercluster can be split)
    int total_vars = 0;
    for (int i = 0; i < CALLSETS*HAPS; i++) {
        total_vars += var_end_indices[i] - var_start_indices[i];
    }
    if (print) printf("total vars: %d\n", total_vars);
    if (total_vars < 2) return var_best_split_indices; // empty

    // get position and hap of next variant
    var_info curr_var = get_next_variant_info(vars, split_indices, var_end_indices);
    if (print) printf("curr var: on hap %d, range %d-%d\n", curr_var.hap_idx, curr_var.start_pos, curr_var.end_pos);
    split_indices[curr_var.hap_idx]++;
    var_info next_var = get_next_variant_info(vars, split_indices, var_end_indices);
    if (print) printf("next var: on hap %d, range %d-%d\n", next_var.hap_idx, next_var.start_pos, next_var.end_pos);
    while (next_var.hap_idx >= 0) {

        // calculate max split size reduction factor
        int gap = std::max(0, next_var.start_pos - curr_var.end_pos);
        double size_reduction_factor = std::max(double((curr_var.end_pos + gap/2) - orig_sc_beg_pos) / orig_sc_size,
                double(orig_sc_end_pos - (curr_var.end_pos + gap/2)) / orig_sc_size);
        double splits_to_halve_size = -1 / log2(size_reduction_factor);

        // calculate overlap (could weight by number of haps overlapping)
        double split_score = gap / splits_to_halve_size;
        if (print) printf("indices: [%d, %d, %d, %d], gap: %d, frac: %f, splits: %f, score: %f\n",
                split_indices[0], split_indices[1], split_indices[2], split_indices[3], gap, size_reduction_factor, splits_to_halve_size, split_score);
        if (split_score > best_split_score) {
            best_split_score = split_score;
            var_best_split_indices = split_indices;
        }

        curr_var = next_var;
        if (print) printf("curr var: on hap %d, range %d-%d\n", curr_var.hap_idx, curr_var.start_pos, curr_var.end_pos);
        split_indices[curr_var.hap_idx]++;
        next_var = get_next_variant_info(vars, split_indices, var_end_indices);
        if (print) printf("next var: on hap %d, range %d-%d\n", next_var.hap_idx, next_var.start_pos, next_var.end_pos);
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
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%sClustering %s VCF%s '%s'", 
            COLOR_PURPLE, callset_strs[callset].data(), 
            COLOR_WHITE, vcf->filename.data());
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
            vars->left_reaches = left_reach;
            vars->right_reaches = right_reach;

            if (vars->clusters[vars->clusters.size()-1] != vars->n) {
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
                int end = std::min(vcf->lengths[ctg_idx], vars->poss[end_idx-1] + vars->rlens[end_idx-1] + 1);
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
					if (beg_pos == 0) break;
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
					if (end_pos == vcf->lengths[ctg_idx]) break;

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
    vars->left_reaches = left_reach;
    vars->right_reaches = right_reach;
}
