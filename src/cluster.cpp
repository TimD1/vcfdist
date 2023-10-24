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
    this->phase.push_back(PHASE_NONE);
    this->orig_phase_dist.push_back(-1);
    this->swap_phase_dist.push_back(-1);
    this->n++;
}

void ctgSuperclusters::set_phase(
           int sc_idx, int phase, 
           int orig_phase_dist, 
           int swap_phase_dist) {
    this->phase[sc_idx] = phase;
    this->orig_phase_dist[sc_idx] = orig_phase_dist;
    this->swap_phase_dist[sc_idx] = swap_phase_dist;
}

/******************************************************************************/

std::vector< std::vector< std::vector<int> > > 
sort_superclusters(std::shared_ptr<superclusterData> sc_data) {

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
                auto vars = ctg_scs->ctg_variants[QUERY][hap];
                if (vars->clusters.size() == 0) continue;
                auto scs = ctg_scs->superclusters[QUERY][hap];
                /* printf("nscs=%d, sc_size=%d, sc_idx=%d\n", ctg_scs->n, int(scs.size()), sc_idx); */
                /* printf("clust_beg=%d, clust_end=%d, clust_size=%d\n", scs[sc_idx], scs[sc_idx+1], int(vars->clusters.size())); */
                int var_beg = vars->clusters[scs[sc_idx]];
                int var_end = vars->clusters[scs[sc_idx+1]];
                /* printf("var_beg=%d, var_end=%d, nvar=%d, var_size=%d\n", var_beg, var_end, vars->n, int(vars->poss.size())); */
                // calculate query len after applying variants
                int query_len = ctg_scs->ends[sc_idx] - ctg_scs->begs[sc_idx];
                for (int var = var_beg; var < var_end; var++) {
                    query_len += (vars->alts[var].size() - vars->refs[var].size());
                }
                max_query_len = std::max(max_query_len, query_len);
            }
            
            // get max truth length
            int max_truth_len = 0;
            for (int hap = 0; hap < HAPS; hap++) {
                // get range of included variants
                auto vars = ctg_scs->ctg_variants[TRUTH][hap];
                if (vars->clusters.size() == 0) continue;
                auto scs = ctg_scs->superclusters[TRUTH][hap];
                int var_beg = vars->clusters[scs[sc_idx]];
                int var_end = vars->clusters[scs[sc_idx+1]];
                // calculate truth len after applying variants
                int truth_len = ctg_scs->ends[sc_idx] - ctg_scs->begs[sc_idx];
                for (int var = var_beg; var < var_end; var++) {
                    truth_len += (vars->alts[var].size() - vars->refs[var].size());
                }
                max_truth_len = std::max(max_truth_len, truth_len);
            }

            // calculate memory usage
            // 10 comes from:
            // (4) aln_ptrs    = 1 byte * 2 haps * 2 phasings
            // (2) path_ptrs   = 1 byte * 2 haps
            // (4) path_scores = 2 byte * 2 haps
            // 2 is a fudge factor that I'm adding for now
            size_t mem = size_t(max_query_len) * size_t(max_truth_len) * 10 * 2;
            double mem_gb = mem / (1000.0 * 1000.0 * 1000.0);
            if (mem_gb > g.max_ram) {
                WARN("Max (%.3fGB) RAM exceeded (%.3fGB req) for supercluster %s:%d-%d, running anyways", 
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

    // generate variant clusters
    if (g.simple_cluster) {
        this->gap_supercluster();
    } else {
        this->supercluster();
    }

    this->transfer_phase_sets();
}


/******************************************************************************/

/* Using per-variant phase set information from PS tags, calculate per-supercluster
 * phase sets. */
void superclusterData::transfer_phase_sets() {
    bool print = false;

    int total_query_phase_sets = 1;
    int total_truth_phase_sets = 1;
    int total_multiple_phase_sets = 0;
    int total_non_increasing = 0;
    std::vector<int> query_phase_set_sizes;
    std::vector<int> truth_phase_set_sizes;
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
                    first_pos = q1v->poss[q1i]; phase_set = q1v->phase_sets[q1i];
                }
                break;
            }
        }
        for (; q2i < q2v->n; q2i++) {
            if (q2v->phase_sets[q2i] != 0) {
                if (q2v->poss[q2i] < first_pos) {
                    first_pos = q2v->poss[q2i]; phase_set = q2v->phase_sets[q2i];
                }
                break;
            }
        }
        for (; t1i < t1v->n; t1i++) {
            if (t1v->phase_sets[t1i] != 0) {
                if (t1v->poss[t1i] < first_pos) {
                    first_pos = t1v->poss[t1i]; phase_set = t1v->phase_sets[t1i];
                }
                break;
            }
        }
        for (; t2i < t2v->n; t2i++) {
            if (t2v->phase_sets[t2i] != 0) {
                if (t2v->poss[t2i] < first_pos) {
                    first_pos = t2v->poss[t2i]; phase_set = t2v->phase_sets[t2i];
                }
                break;
            }
        }
        if (print) printf("first phase set: %d\n", phase_set);

        // carry phase_set across superclusters (if no heterozygous variants)
        int query_phase_set = phase_set;
        int truth_phase_set = phase_set;
        int query_ps_beg = phase_set; int query_ps_end = 0;
        int truth_ps_beg = phase_set; int truth_ps_end = 0;
        for (int sci = 0; sci < ctg_scs->n; sci++) { // for each supercluster
            int query_ps_ct = 0;
            int truth_ps_ct = 0;
            if (print) printf("supercluster: %d", sci);

            // check all variants in supercluster for each haplotype
            if (q1v->n)
            for (q1i = q1v->clusters[q1sc[sci]]; // QUERY HAP 1
                    q1i < q1v->clusters[q1sc[sci+1]]; q1i++) {
                if (q1v->phase_sets[q1i]) { // non-zero, has PS tag
                    if (q1v->phase_sets[q1i] > query_phase_set) { // new PS
                        phase_set = query_phase_set = q1v->phase_sets[q1i];
                        total_query_phase_sets++; query_ps_ct++;
                        query_phase_set_sizes.push_back(query_ps_end - query_ps_beg); // add phase set
                        if (print) printf(" Q1 PS:%d %d-%d", phase_set, query_ps_beg, query_ps_end);
                        query_ps_beg = q1v->poss[q1i];
                        query_ps_end = q1v->poss[q1i] + q1v->rlens[q1i];
                    } else if (q1v->phase_sets[q1i] == query_phase_set) { // same 
                        if (!query_ps_ct) { query_ps_ct++; }
                        query_ps_end = std::max(query_ps_end, q1v->poss[q1i] + q1v->rlens[q1i]);
                    } else {
                        total_non_increasing++;
                    }
                }
            }
            if (q2v->n)
            for (q2i = q2v->clusters[q2sc[sci]]; // QUERY HAP 2
                    q2i < q2v->clusters[q2sc[sci+1]]; q2i++) {
                if (q2v->phase_sets[q2i]) { // non-zero, has PS tag
                    if (q2v->phase_sets[q2i] > query_phase_set) { // new
                        phase_set = query_phase_set = q2v->phase_sets[q2i];
                        total_query_phase_sets++; query_ps_ct++;
                        query_phase_set_sizes.push_back(query_ps_end - query_ps_beg); // add phase set
                        if (print) printf(" Q2 PS:%d %d-%d", phase_set, query_ps_beg, query_ps_end);
                        query_ps_beg = q2v->poss[q2i];
                        query_ps_end = q2v->poss[q2i] + q2v->rlens[q2i];
                    } else if (q2v->phase_sets[q2i] == query_phase_set) { // same 
                        if (!query_ps_ct) { query_ps_ct++; }
                        query_ps_end = std::max(query_ps_end, q2v->poss[q2i] + q2v->rlens[q2i]);
                    } else {
                        total_non_increasing++;
                    }
                }
            }
            if (t1v->n)
            for (t1i = t1v->clusters[t1sc[sci]]; // TRUTH HAP 1
                    t1i < t1v->clusters[t1sc[sci+1]]; t1i++) {
                if (t1v->phase_sets[t1i]) { // non-zero, has PS tag
                    if (t1v->phase_sets[t1i] > truth_phase_set) { // new
                        phase_set = truth_phase_set = t1v->phase_sets[t1i];
                        total_truth_phase_sets++; truth_ps_ct++;
                        truth_phase_set_sizes.push_back(truth_ps_end - truth_ps_beg); // add phase set
                        if (print) printf(" T1 PS:%d %d-%d", phase_set, truth_ps_beg, truth_ps_end);
                        truth_ps_beg = t1v->poss[t1i];
                        truth_ps_end = t1v->poss[t1i] + t1v->rlens[t1i];
                    } else if (t1v->phase_sets[t1i] == truth_phase_set) { // same 
                        if (!truth_ps_ct) { truth_ps_ct++; }
                        truth_ps_end = std::max(truth_ps_end, t1v->poss[t1i] + t1v->rlens[t1i]);
                    } else {
                        total_non_increasing++;
                    }
                }
            }
            if (t2v->n)
            for (t2i = t2v->clusters[t2sc[sci]]; // TRUTH HAP 2
                    t2i < t2v->clusters[t2sc[sci+1]]; t2i++) {
                if (t2v->phase_sets[t2i]) { // non-zero, has PS tag
                    if (t2v->phase_sets[t2i] > truth_phase_set) { // new
                        phase_set = truth_phase_set = t2v->phase_sets[t2i];
                        total_truth_phase_sets++; truth_ps_ct++;
                        truth_phase_set_sizes.push_back(truth_ps_end - truth_ps_beg); // add phase set
                        if (print) printf(" T2 PS:%d %d-%d", phase_set, truth_ps_beg, truth_ps_end);
                        truth_ps_beg = t2v->poss[t2i];
                        truth_ps_end = t2v->poss[t2i] + t2v->rlens[t2i];
                    } else if (t2v->phase_sets[t2i] == truth_phase_set) { // same 
                        if (!truth_ps_ct) { truth_ps_ct++; }
                        truth_ps_end = std::max(truth_ps_end, t2v->poss[t2i] + t2v->rlens[t2i]);
                    } else {
                        total_non_increasing++;
                    }
                }
            }

            // if one supercluster contains variants in multiple phase sets
            if (query_ps_ct > 1 || truth_ps_ct > 1) {
                total_multiple_phase_sets++;
            }

            // add phase_set to supercluster
            ctg_scs->phase_sets.push_back(phase_set);
            if (print) printf(" final=%d\n", phase_set);
        }

        // add final phase set on contig
        query_phase_set_sizes.push_back(query_ps_end - query_ps_beg);
        truth_phase_set_sizes.push_back(truth_ps_end - truth_ps_beg);
    }

    // calculate phaseset NG50
    size_t total_bases = 0;
    for (size_t i = 0; i < this->contigs.size(); i++) {
        total_bases += lengths[i];
    }

    int query_ps_ng50 = calc_ng50(query_phase_set_sizes, total_bases);
    int truth_ps_ng50 = calc_ng50(truth_phase_set_sizes, total_bases);

    // print
    if (total_multiple_phase_sets)
        WARN("%d total superclusters contain variants from multiple phase sets (PS)",
                total_multiple_phase_sets)

    if (total_non_increasing)
        WARN("%d total phase sets (PS) are in non-increasing order",
                total_non_increasing)

    if (g.verbosity >= 1) INFO("        Input QUERY phase sets: %d", total_query_phase_sets);
    if (g.verbosity >= 1) INFO("     Input QUERY phaseset NG50: %d", query_ps_ng50);
    if (g.verbosity >= 1) INFO("        Input TRUTH phase sets: %d", total_truth_phase_sets);
    if (g.verbosity >= 1) INFO("     Input TRUTH phaseset NG50: %d", truth_ps_ng50);
}


/* Supercluster using left_reach and right_reach of each cluster, calculated
 * during wf_swg_cluster().
 */
void superclusterData::supercluster() {
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
        auto vars = this->superclusters[ctg]->ctg_variants;
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
                if (vars[i>>1][i&1]->clusters.size())
                    this_vars += vars[i>>1][i&1]->clusters[ next_brks[i] ] - 
                        vars[i>>1][i&1]->clusters[ brks[i] ];
            }
            most_vars = std::max(most_vars, this_vars);
            total_vars += this_vars;

            // debug print
            if (false) {
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
                printf("variants: %d, bases: %d\n", this_vars, end_pos-beg_pos);
            }

            // save alignment information
            this->superclusters[ctg]->add_supercluster(brks, beg_pos, end_pos);

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


/* Supercluster by enforcing minimum gap between clusters.
 */
void superclusterData::gap_supercluster() {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[4/8] Gap superclustering TRUTH and QUERY variants%s",
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
        for (int i = 0; i < CALLSETS*HAPS; i++)
            nvars += this->superclusters[ctg]->ctg_variants[i>>1][i&1]->n;
        if (!nvars) continue;

        // for each cluster of variants (merge query and truth haps)
        auto vars = this->superclusters[ctg]->ctg_variants;
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
                this_vars += vars[i>>1][i&1]->clusters[ next_brks[i] ] - 
                    vars[i>>1][i&1]->clusters[ brks[i] ];
            }
            most_vars = std::max(most_vars, this_vars);
            total_vars += this_vars;

            // debug print
            if (false) {
                printf("\nSUPERCLUSTER: %d\n", this->superclusters[ctg]->n);
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
            this->superclusters[ctg]->add_supercluster(brks, beg_pos, end_pos);

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

/* Cluster variants based on minimum gap length for independence */
void gap_cluster(std::shared_ptr<variantData> vcf, int callset) {
    /* Add single-VCF cluster indices to `variantData` */
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%sGap clustering %s VCF%s '%s'", 
            COLOR_PURPLE, callset_strs[callset].data(), 
            COLOR_WHITE, vcf->filename.data());

    // cluster each contig
    int largest_cluster_vars = 0;
    int largest_cluster_bases = 0;
    for (std::string ctg : vcf->contigs) {

        // cluster per-haplotype variants: vcf->variants[hap]
        for (int hap = 0; hap < HAPS; hap++) {
            int prev_end = -g.cluster_min_gap * 2;
            int cluster_start = 0;
            int pos = 0;
            int end = 0;
            int var_idx = 0;
            int prev_var_idx = 0;
            for (; var_idx < vcf->variants[hap][ctg]->n; var_idx++) {
                pos = vcf->variants[hap][ctg]->poss[var_idx];
                end = pos + vcf->variants[hap][ctg]->rlens[var_idx];
                if (pos - prev_end > g.cluster_min_gap) {
                    // update summary metrics
                    largest_cluster_bases = std::max(largest_cluster_bases, prev_end-cluster_start);
                    largest_cluster_vars = std::max(largest_cluster_vars, var_idx - prev_var_idx);
                    prev_var_idx = var_idx;
                    cluster_start = pos;
                    // add cluster
                    vcf->variants[hap][ctg]->add_cluster(var_idx);
                }
                prev_end = std::max(prev_end, end);
            }
            vcf->variants[hap][ctg]->add_cluster(var_idx);
        }
    }
    if (g.verbosity >= 1) INFO("  Largest cluster (bases): %d", largest_cluster_bases);
    if (g.verbosity >= 1) INFO("  Largest cluster  (vars): %d", largest_cluster_vars);
}


/******************************************************************************/


/* Add single-VCF cluster indices to `variantData`. This version assumes that
 * all variant calls are true positives (doesn't allow skipping)
 */
void wf_swg_cluster(variantData * vcf, std::string ctg, int hap,
        int sub, int open, int extend) {
    bool print = false;

    // allocate this memory once, use on each cluster
    std::vector<int> offs_buffer(MATS * g.max_size*2 * 
            std::max(sub, open+extend), -2);

    // init: each variant is its own cluster
    auto vars = vcf->variants[hap][ctg];
    int nvar = vars->n;
    if (!nvar) return;

    // mark variant boundary between clusters (hence n+1), 
    // and whether it was active on the previous iteration.
    // initially, one variant per cluster
    std::vector<int> prev_clusters(nvar+1);
    for (int i = 0; i < nvar+1; i++) 
        prev_clusters[i] = i;
    std::vector<bool> prev_active(nvar+1, true);

    std::vector<int> right_reach(nvar+1), left_reach(nvar+1);
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

        // set far left/right to avoid OOB for first/last clusters
        left_reach[0] = vars->poss[0];
        right_reach[prev_clusters.size()-2] = vars->poss[nvar-1];
        // set sentinel reaches
        left_reach[prev_clusters.size()-1] = std::numeric_limits<int>::max();
        right_reach[prev_clusters.size()-1] = std::numeric_limits<int>::max();

        // update all cluster reaches
        for (size_t clust = 0; clust < prev_clusters.size(); clust++) {

            // only compute necessary reaches (adjacent merge)
            bool left_compute = prev_active[clust];
            bool right_compute = prev_active[clust];

            // no left cluster, or sentinel
            if (clust == 0 || clust == prev_clusters.size()-1) {
                left_compute = false;
            }
            // no right cluster, or sentinel
            if (clust >= prev_clusters.size()-2) {
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
                while (reach == ref_len-1) { // iterative doubling
                    ref_len *= 2;
                    int beg_pos = end_pos - std::max(0, main_diag) - ref_len - score/extend-3;
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
                while (reach == ref_len-1) { // iterative doubling
                    ref_len *= 2;
                    int end_pos = beg_pos + std::max(0, main_diag) + ref_len + score/extend+3;
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
        /* tmp_left_reach.push_back(std::numeric_limits<int>::max()); */
        /* tmp_right_reach.push_back(std::numeric_limits<int>::min()); */
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
