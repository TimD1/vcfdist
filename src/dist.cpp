#include <cassert>
#include <chrono>
#include <cstdio>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "cluster.h"
#include "dist.h"
#include "print.h"

template <typename T>
inline bool contains(const std::unordered_set<T> & wave, const T & idx) {
    return wave.find(idx) != wave.end();
}
template <typename T, typename U>
inline bool contains(const std::unordered_map<T,U> & wave, const T & idx) {
    return wave.find(idx) != wave.end();
}

/******************************************************************************/


int calc_ng50(std::vector<int> phase_blocks, size_t total_bases) {
    size_t total_phased = 0;
    std::sort(phase_blocks.begin(), phase_blocks.end(), std::greater<>());
    for (size_t i = 0; i < phase_blocks.size(); i++) {
        total_phased += phase_blocks[i];
        if (total_phased >= total_bases / 2)
            return phase_blocks[i];
    }
    return 0;
}


/******************************************************************************/


/* Generate the new sequence by applying variants to the reference. */
std::string generate_str(
        std::shared_ptr<fastaData> ref, 
        std::shared_ptr<ctgVariants> vars, const std::string & ctg,
        int beg_idx, int end_idx, int beg_pos, int end_pos, 
        int min_qual /* = 0 */) {

    int var_idx = beg_idx;
    std::string str = "";
    while (var_idx < vars->n && vars->poss[var_idx] < beg_pos) var_idx++;
    for (int ref_pos = beg_pos; ref_pos < end_pos; ) {

        // VARIANT
        if (var_idx < end_idx && ref_pos == vars->poss[var_idx]) {
            if (vars->var_quals[var_idx] >= min_qual) {
                switch (vars->types[var_idx]) {
                    case TYPE_INS:
                        str += vars->alts[var_idx];
                        break;
                    case TYPE_DEL:
                        ref_pos += vars->refs[var_idx].size();
                        break;
                    case TYPE_SUB:
                        str += vars->alts[var_idx];
                        ref_pos++;
                        break;
                    case TYPE_CPX:
                        str += vars->alts[var_idx];
                        ref_pos += vars->refs[var_idx].size();
                        break;
                }
            }
            var_idx++; // next variant

        } else { // NO VARIANT
            try {
                // add entire ref seq up to next variant or region end
                int ref_end = (var_idx < end_idx) ? 
                    std::min(end_pos, vars->poss[var_idx]) : end_pos;
                if (ref_end < ref_pos) {
                    printf("ref_pos:%d ref_end:%d beg_pos: %d end_pos: %d\n",
                            ref_pos, ref_end, beg_pos, end_pos);
                    for(int i = beg_idx; i < end_idx; i++) {
                        printf("%s:%d %s -> %s\n", ctg.data(), vars->poss[i], 
                                vars->refs[i].data(), vars->alts[i].data());
                    }
                    ERROR("No variant, but ref_end < ref_pos (generate_str)");
                }
                str += ref->fasta.at(ctg).substr(ref_pos, ref_end-ref_pos);
                ref_pos = ref_end;
            } catch (const std::out_of_range & e) {
                ERROR("Contig '%s' not in reference FASTA or position out of range (generate_str)", ctg.data());
            }
        }
    }
    return str;
}


/******************************************************************************/


/* Calculate initial forward-pass alignment for truth string and query graph.
 * This function generates a pointer matrix and alignment score.
 */
int calc_prec_recall_aln(
        const std::shared_ptr<Graph> graph,
        std::unordered_map<idx4, idx4> & ptrs,
        bool print
        ) {

    // init matrices
    std::unordered_set<idx4> done; // visited cells
    std::queue<idx4> queue; // still to be explored in this wave
    
    // set first wavefront
    idx4 start(0, 0, 0, 0);
    queue.push(start);
    done.insert(start);
    ptrs[start] = idx4(0, 0, -1, -1);
    int score = 0;

    // continue looping until full alignment found
    std::unordered_set<idx4> curr_wave; // everything explored this wave
    std::unordered_set<idx4> prev_wave; // everything explored prev wave
    while (score <= g.max_dist) {
        /* if (print) printf("  score = %d\n", score); */
        if (queue.empty()) ERROR("Empty queue in 'prec_recall_aln()'.");

        // EXTEND WAVEFRONT (stay at same score)
        while (!queue.empty()) {
            idx4 x = queue.front(); queue.pop();
            /* if (print) printf("    x = node (%d, %d) cell (%d, %d)\n", x.qni, x.tni, x.qi, x.ti); */
            prev_wave.insert(x);

            // allow match on query
            idx4 y(x.qni, x.tni, x.qi, x.ti);
            bool match = false;
            while (y.qi+1 < int(graph->qseqs[y.qni].length()) &&
                   y.ti+1 < int(graph->tseqs[y.tni].length()) &&
                   graph->qseqs[y.qni][y.qi+1] == graph->tseqs[y.tni][y.ti+1]) {
                y.qi++;
                y.ti++;
                match = true;
            }
            if (match && !contains(done, y) && !contains(curr_wave, y)) {
                /* if (print) printf("      y = node (%d, %d) cell (%d, %d)\n", y.qni, y.tni, y.qi, y.ti); */
                queue.push(y); curr_wave.insert(y); ptrs[y] = x;
            }

            // allow bottom-right corner to move diagonally into next truth and query nodes
            // NOTE: separating this case out allows sync points between adjacent variants
            if (x.qi == int(graph->qseqs[x.qni].length())-1 && x.ti == int(graph->tseqs[x.tni].length())-1) {
                for (int qni : graph->qnexts[x.qni]) { // for all next nodes
                    for (int tni : graph->tnexts[x.tni]) {
                        idx4 z(qni, tni, 0, 0);
                        if (!contains(done, z) && !contains(curr_wave, z)) {
                            /* if (print) printf("      z = node (%d, %d) cell (%d, %d)\n", z.qni, z.tni, z.qi, z.ti); */
                            queue.push(z); curr_wave.insert(z); ptrs[z] = x;
                        }
                    }
                }
            }

            // allow last row to move into first row of all next query nodes
            if (x.qi == int(graph->qseqs[x.qni].length())-1 && x.ti < int(graph->tseqs[x.tni].length())) {
                for (int qni : graph->qnexts[x.qni]) { // for all next nodes
                    idx4 z(qni, x.tni, 0, x.ti);
                    if (!contains(done, z) && !contains(curr_wave, z)) {
                        /* if (print) printf("      z = node (%d, %d) cell (%d, %d)\n", z.qni, z.tni, z.qi, z.ti); */
                        queue.push(z); curr_wave.insert(z); ptrs[z] = x;
                    }
                }
            }

            // allow last col to move into first col of next truth node
            if (x.ti == int(graph->tseqs[x.tni].length())-1 && x.qi < int(graph->qseqs[x.qni].length())) {
                for (int tni : graph->tnexts[x.tni]) { // for all next nodes
                    idx4 z(x.qni, tni, x.qi, 0);
                    if (!contains(done, z) && !contains(curr_wave, z)) {
                        /* if (print) printf("      z = node (%d, %d) cell (%d, %d)\n", z.qni, z.tni, z.qi, z.ti); */
                        queue.push(z); curr_wave.insert(z); ptrs[z] = x;
                    }
                }
            }
        }

        // mark all cells visited this wave as done
        for (idx4 x : curr_wave) { done.insert(x); }
        curr_wave.clear();

        // exit if we're done aligning
        idx4 end = idx4(graph->qnodes-1, graph->tnodes-1, 
                    int(graph->qseqs[graph->qnodes-1].length()-1), 
                    int(graph->tseqs[graph->tnodes-1].length()-1));
        if (contains(done, end)) break;

        // NEXT WAVEFRONT (increase score by one)
        for (idx4 x : prev_wave) {
            if (x.qi+1 < int(graph->qseqs[x.qni].length())) { // INS
                idx4 y(x.qni, x.tni, x.qi+1, x.ti);
                if (!contains(done, y) && !contains(curr_wave, y)) {
                    queue.push(y); curr_wave.insert(y); ptrs[y] = x;
                }
            }
            if (x.ti+1 < int(graph->tseqs[x.tni].length())) { // DEL
                idx4 y(x.qni, x.tni, x.qi, x.ti+1);
                if (!contains(done, y) && !contains(curr_wave, y)) {
                    queue.push(y); curr_wave.insert(y); ptrs[y] = x;
                }
            }
            if (x.qi+1 < int(graph->qseqs[x.qni].length()) && 
                    x.ti+1 < int(graph->tseqs[x.tni].length())) { // SUB
                idx4 y(x.qni, x.tni, x.qi+1, x.ti+1);
                if (!contains(done, y) && !contains(curr_wave, y)) {
                    queue.push(y); curr_wave.insert(y); ptrs[y] = x;
                }
            }
        }
        prev_wave.clear();
        score++;
    }

    if (print) {
        printf("Alignment score: %d\n", score);
        printf("Alignment:\n");
        print_graph_ptrs(graph, ptrs);
    }

    return score;
}


/******************************************************************************/


void evaluate_variants(std::shared_ptr<ctgSuperclusters> scs, int sc_idx,
			std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hi, bool print) {

    bool done = false;
    while (not done) {

        // graph is constructed only from unevaluated variants
        std::shared_ptr<Graph> graph(new Graph(scs, sc_idx, ref, ctg, truth_hi));
        std::shared_ptr<ctgVariants> qvars = graph->sc->callset_vars[QUERY];
        std::shared_ptr<ctgVariants> tvars = graph->sc->callset_vars[TRUTH];
        if (print) graph->print();

        std::unordered_map<idx4, idx4> ptrs;
        int score = calc_prec_recall_aln(graph, ptrs, print);
        bool aligned = score <= g.max_dist;
        if (aligned) { // alignment succeeded
            calc_prec_recall(graph, ptrs, truth_hi, print);
            done = true;
        }
        ptrs.clear();

        // NOTE: if alignment failed because it was too expensive, this can only be caused by a FN
        // truth variant, since query variants can be skipped and the reference sections match.
        // Try alignments without some of the largest truth variants

        // NOTE: if we did align successfully and found a FN INDEL, re-evaluate all FP and FN variants 
        // besides the largest FN. The motivation for this is that large FN truth SVs will have a sync 
        // group that extends pretty far left and right, "swallowing" other correct variant calls. 
        // Since correctness is determined per sync group, many TP SNP calls can't be identified unless
        // evaluated without the presence of this FN SV. So we remove it and re-evaluate.
        std::vector< std::pair<int, int> > exclude_sizes; // pairs of (truth_var_size, truth_var_idx)
        for (int tni = 0; tni < graph->tnodes; tni++) {
            if (graph->ttypes[tni] != TYPE_REF) {
                int tvar_idx = graph->tidxs[tni];
                if (not aligned or (tvars->errtypes[truth_hi][tvar_idx] == ERRTYPE_FN
                            and tvars->types[tvar_idx] != TYPE_SUB)) {
                    done = false;
                    int tvar_size = std::max(tvars->alts[tvar_idx].size(), tvars->refs[tvar_idx].size());
                    exclude_sizes.push_back(std::make_pair(tvar_size, tvar_idx));
                }
            }
        }
        // NOTE: sort FNs (if aligned) or all truth vars (if not aligned) by decreasing variant size
        std::sort(exclude_sizes.rbegin(), exclude_sizes.rend());
        if (print) for (int i = 0; i < int(exclude_sizes.size()); i++) {
            int tvar_idx = exclude_sizes[i].second;
            printf("  exclude size %d, var %d = %s:%d (%s,%s)\n",
                    exclude_sizes[i].first, tvar_idx, ctg.data(), 
                    tvars->poss[tvar_idx], tvars->refs[tvar_idx].data(),
                    tvars->alts[tvar_idx].data());
        }
        // try a number of alignments without some of the largest variants, choose the best one
        if (not done) { // there are potential variants to exclude
            std::vector< std::pair<int, int> > exclude_dists;
            for (int retry = 0; retry < std::min(g.max_retries, int(exclude_sizes.size())); retry++) {

                // retry with all variants except the chosen variant to exclude (FN)
                for (int tni = 0; tni < graph->tnodes; tni++) {
                    if (graph->ttypes[tni] != TYPE_REF) {
                        int tvar_idx = graph->tidxs[tni];
                        if (tvar_idx == exclude_sizes[retry].second) {
                            tvars->errtypes[truth_hi][tvar_idx] = ERRTYPE_FN;
                        } else {
                            tvars->errtypes[truth_hi][tvar_idx] = ERRTYPE_UN;
                        }
                    }
                }
                for (int qni = 0; qni < graph->qnodes; qni++) {
                    if (graph->qtypes[qni] != TYPE_REF) {
                        int qvar_idx = graph->qidxs[qni];
                        qvars->errtypes[truth_hi][qvar_idx] = ERRTYPE_UN;
                    }
                }

                // retry alignment (with one large variant excluded, as well as all TPs)
                std::shared_ptr<Graph> retry_graph(new Graph(scs, sc_idx, ref, ctg, truth_hi));
                if (print) retry_graph->print();
                std::unordered_map<idx4, idx4> retry_ptrs;
                int dist = calc_prec_recall_aln(retry_graph, retry_ptrs, print);
                exclude_dists.push_back(std::make_pair(dist, exclude_sizes[retry].second));
            }

            // sort retried alignments in order of resulting edit distance
            std::sort(exclude_dists.begin(), exclude_dists.end());
            if (print) for (int i = 0; i < int(exclude_dists.size()); i++) {
                int tvar_idx = exclude_dists[i].second;
                printf("  exclude dist %d, var %d = %s:%d (%s,%s)\n",
                        exclude_dists[i].first, tvar_idx, ctg.data(), 
                        tvars->poss[tvar_idx], tvars->refs[tvar_idx].data(),
                        tvars->alts[tvar_idx].data());
            }

            // the variant for which excluding resulted in lowest edit dist should be considered FN
            // prepare variants for the next iteration
            for (int tni = 0; tni < graph->tnodes; tni++) {
                if (graph->ttypes[tni] != TYPE_REF) {
                    int tvar_idx = graph->tidxs[tni];
                    if (tvar_idx == exclude_dists[0].second) { // remove this variant from eval
                        tvars->errtypes[truth_hi][tvar_idx] = ERRTYPE_FN;
                    } else {
                        tvars->errtypes[truth_hi][tvar_idx] = ERRTYPE_UN;
                        tvars->sync_group[truth_hi][tvar_idx] = 0;
                        tvars->callq[truth_hi][tvar_idx] = 0;
                        tvars->ref_ed[truth_hi][tvar_idx] = 0;
                        tvars->query_ed[truth_hi][tvar_idx] = 0;
                        tvars->credit[truth_hi][tvar_idx] = 0;
                    }
                }
            }
            for (int qni = 0; qni < graph->qnodes; qni++) {
                if (graph->qtypes[qni] != TYPE_REF) {
                    int qvar_idx = graph->qidxs[qni];
                    qvars->set_var_calcgt_on_hap(qvar_idx, truth_hi, false, true);
                    qvars->errtypes[truth_hi][qvar_idx] = ERRTYPE_UN;
                    qvars->sync_group[truth_hi][qvar_idx] = 0;
                    qvars->callq[truth_hi][qvar_idx] = 0;
                    qvars->ref_ed[truth_hi][qvar_idx] = 0;
                    qvars->query_ed[truth_hi][qvar_idx] = 0;
                    qvars->credit[truth_hi][qvar_idx] = 0;
                }
            }
        }
    }
}


/******************************************************************************/


void calc_prec_recall(
        const std::shared_ptr<Graph> graph,
        const std::unordered_map<idx4, idx4> & ptrs,
        int truth_hap, bool print
        ) {

    idx4 end(graph->qnodes-1, graph->tnodes-1,
            graph->qseqs[graph->qnodes-1].length()-1,
            graph->tseqs[graph->tnodes-1].length()-1);
    idx4 curr = end;
    int prev_query_ref_pos = graph->ref.size();
    int prev_truth_pos = graph->truth.size();
    int sync_group = 0;
    int query_dist = 0;
    std::vector<int> sync_tvars;
    std::vector<int> sync_qvars;
    std::shared_ptr<ctgVariants> qvars = graph->sc->callset_vars[QUERY];
    std::shared_ptr<ctgVariants> tvars = graph->sc->callset_vars[TRUTH];

    // all query variants default to FP unless path chosen
    for (int qni = 0; qni < graph->qnodes; qni++) {
        if (graph->qtypes[qni] != TYPE_REF) {
            int qvar_idx = graph->qidxs[qni];
            qvars->errtypes[truth_hap][qvar_idx] = ERRTYPE_FP;
            qvars->callq[truth_hap][qvar_idx] = qvars->var_quals[qvar_idx];
        }
    }

    while (curr != idx4(0, 0, -1, -1)) {

        idx4 prev = ptrs.at(curr);

        if (print) printf("curr = node (%d, %d) cell (%d, %d) poss (%d, %d)",
                curr.qni, curr.tni, curr.qi, curr.ti, 
                graph->qbegs[curr.qni]+curr.qi, graph->tbegs[curr.tni]+curr.ti);
        if (print) printf("\tprev = node (%d, %d) cell (%d, %d)\n",
                prev.qni, prev.tni, prev.qi, prev.ti);

        /* printf("%s submatrix, %s diag, query %s, truth %s, %s pos\n", */
        /*     (prev.qni == curr.qni && prev.tni == curr.tni) ? "same" : "diff", */
        /*     (prev.qi+1 == curr.qi && prev.ti+1 == curr.ti) ? "is" : "is not", */
        /*     (graph->qtypes[curr.qni] == TYPE_REF) ? "ref" : "var", */
        /*     (graph->ttypes[curr.tni] == TYPE_REF) ? "ref" : "var", */
        /*     (graph->tbegs[curr.tni] + curr.ti == graph->qbegs[curr.qni] + curr.qi) ? "same" : "diff"); */

        // if we move out of a new query variant, mark it as included
        if (prev.qni != curr.qni && // new query node
                graph->qtypes[curr.qni] != TYPE_REF) { // node is variant
            int qvar_idx = graph->qidxs[curr.qni];
            sync_qvars.push_back(qvar_idx);
            if (print) printf("new query variant: %d\n", qvar_idx);
        }
        // if we move out of a query variant, include it in sync group and ref dist calc
        if (prev.tni != curr.tni && // new truth node
                graph->ttypes[curr.tni] != TYPE_REF) { // node is variant
            int tvar_idx = graph->tidxs[curr.tni];
            sync_tvars.push_back(tvar_idx);
            if (print) printf("new truth variant: %d\n", tvar_idx);
        } 
        // if the alignment is a substitution, insertion, or deletion
        if (prev.qni == curr.qni && prev.tni == curr.tni && // same matrix
                (prev.qi == curr.qi || prev.ti == curr.ti || // insertion or deletion
                 graph->tseqs[curr.tni][curr.ti] != graph->qseqs[curr.qni][curr.qi]) // substitution
                ) {
            if (print) printf("non-match step\n");
            query_dist += 1;
        }

        // check if this movement is a sync point (ref main diag mvmt or between var main diag mvmt)
        int truth_ref_pos = graph->tbegs[curr.tni] + curr.ti;
        int query_ref_pos = graph->qbegs[curr.qni] + curr.qi;
        bool on_main_diag = (truth_ref_pos == query_ref_pos);
        bool same_submatrix = curr.qni == prev.qni && curr.tni == prev.tni;
        bool diff_submatrix = curr.qni != prev.qni && curr.tni != prev.tni; // both must differ
        bool ref_query_move = graph->qtypes[curr.qni] == TYPE_REF 
            && prev.qni == curr.qni
            && prev.qi < curr.qi;
        bool ref_truth_move = graph->ttypes[curr.tni] == TYPE_REF 
            && prev.tni == curr.tni
            && prev.ti < curr.ti;
        bool sync_point = on_main_diag && (
                (same_submatrix && ref_query_move && ref_truth_move) || diff_submatrix);
        if (sync_point) {
            if (print) printf("potential sync point\n");

            // add sync point
            if (sync_tvars.size() || sync_qvars.size()) {

                // calculate edit distance without variants
                int ref_dist = 0;
                int truth_pos = graph->get_truth_pos(curr.tni, curr.ti);
                wf_ed(graph->ref.substr(query_ref_pos, prev_query_ref_pos - query_ref_pos), 
                        graph->truth.substr(truth_pos, prev_truth_pos - truth_pos),
                        ref_dist);
                if (print) printf("syncing\n");
                float credit = 0;
                if (ref_dist == 0) {
                    // false positive with no nearby truth variant, credit = 0
                } else {
                    credit = 1 - float(query_dist) / ref_dist;
                }
                uint8_t errtype = (credit >= g.credit_threshold) ? ERRTYPE_TP : ERRTYPE_FP;
                float qual = g.max_qual;
                for (int qvar_idx : sync_qvars) {
                    qual = std::min(qual, qvars->var_quals[qvar_idx]);
                }
                for (int qvar_idx : sync_qvars) {
                    if (print) printf("QUERY var %d: BD=%s, BC=%f, Q=%f, SG=%d, RD=%d, QD=%d\n",
                            qvar_idx, error_strs[errtype].data(), credit, 
                            qual, sync_group, ref_dist, query_dist);
                    qvars->errtypes[truth_hap][qvar_idx] = errtype;
                    qvars->callq[truth_hap][qvar_idx] = qual;
                    qvars->sync_group[truth_hap][qvar_idx] = sync_group;
                    qvars->ref_ed[truth_hap][qvar_idx] = ref_dist;
                    qvars->query_ed[truth_hap][qvar_idx] = query_dist;
                    qvars->credit[truth_hap][qvar_idx] = credit;
                    if (errtype == ERRTYPE_TP)
                        qvars->set_var_calcgt_on_hap(qvar_idx, truth_hap, true);
                }
                if (errtype == ERRTYPE_FP) errtype = ERRTYPE_FN;
                for (int tvar_idx : sync_tvars) {
                    if (print) printf("TRUTH var %d: BD=%s, BC=%f, Q=%f, SG=%d, RD=%d, QD=%d\n",
                            tvar_idx, error_strs[errtype].data(), credit, 
                            qual, sync_group, ref_dist, query_dist);
                    tvars->errtypes[truth_hap][tvar_idx] = errtype;
                    tvars->callq[truth_hap][tvar_idx] = qual;
                    tvars->sync_group[truth_hap][tvar_idx] = sync_group;
                    tvars->ref_ed[truth_hap][tvar_idx] = ref_dist;
                    tvars->query_ed[truth_hap][tvar_idx] = query_dist;
                    tvars->credit[truth_hap][tvar_idx] = credit;
                }

                // reset sync group
                sync_group++;
                sync_qvars.clear();
                sync_tvars.clear();
                prev_query_ref_pos = query_ref_pos;
                prev_truth_pos = truth_pos;
                query_dist = 0;
            } 
        }
        curr = prev;
    }
}


/******************************************************************************/


void wf_swg_align(
        const std::string & query, const std::string & truth, 
        int & s, int x, int o, int e, bool print
        ) {

    // init
    int query_len = query.size();
    int truth_len = truth.size();
    int mat_len = query_len + truth_len - 1;
    bool done = false;
    std::vector< std::vector< std::vector<int> > > offs(MATS);
    for (int m = 0; m < MATS; m++) {
        offs[m].push_back(std::vector<int>(mat_len, -2));
    }
    s = 0;
    offs[MAT_SUB][s][query_len-1] = -1; // main diag

    while (true) {

        // EXTEND WAVEFRONT (leave INS, DEL)
        for (int m = MAT_INS; m < MATS; m++) {
            for (int d = 0; d < mat_len; d++) {
                int off = offs[m][s][d];
                int diag = d + 1 - query_len;

                if (off >= 0 && off < query_len &&
                        diag+off >= 0 && diag+off < truth_len &&
                        offs[m][s][d] >= offs[MAT_SUB][s][d]) {
                    offs[MAT_SUB][s][d] = offs[m][s][d];
                    if(print) printf("(S, %d, %d) swap\n", offs[MAT_SUB][s][d], 
                            offs[MAT_SUB][s][d]+d+1-query_len);
                }
            }
        }

        // EXTEND WAVEFRONT (diag, SUB only)
        for (int d = 0; d < mat_len; d++) {
            int off = offs[MAT_SUB][s][d];
            int diag = d + 1 - query_len;

            // extend
            while ( off != -2 && diag + off >= -1 && 
                    off < query_len - 1 && 
                    diag + off < truth_len - 1) {
                if (query[off+1] == truth[diag+off+1]) off++;
                else break;
            }
            if (off > offs[MAT_SUB][s][d])
                if(print) printf("(S, %d, %d) extend\n", off, off+diag);
            offs[MAT_SUB][s][d] = off;

            // finish if done
            if (off == query_len - 1 && off + diag == truth_len - 1)
            { done = true; break; }
        }
        if (done) break;

        // debug print
        for (int mi = 0; mi < MATS; mi++) {
            if(print) printf("\n%s matrix\n", type_strs[mi+1].data());
            if(print) printf("offs %d:", s);
            for (int di = 0; di < int(query.size() + truth.size()-1); di++) {
                if(print) printf("\t%d", offs[mi][s][di]);
            }
            if(print) printf("\n");
        }

        // NEXT WAVEFRONT
        s++;
        if(print) printf("\nscore = %d\n", s);
        for (int m = 0; m < MATS; m++) {
            offs[m].push_back(std::vector<int>(mat_len, -2));
        }

        for (int d = 0; d < mat_len; d++) {
            int diag = d + 1 - query_len;

            // sub (in SUB)
            if (s-x >= 0 && offs[MAT_SUB][s-x][d] != -2 && 
                            offs[MAT_SUB][s-x][d]+1 < query_len &&
                     diag + offs[MAT_SUB][s-x][d]+1 < truth_len &&
                            offs[MAT_SUB][s-x][d]+1 >= offs[MAT_SUB][s][d]) {
                offs[MAT_SUB][s][d] = offs[MAT_SUB][s-x][d] + 1;
                if(print) printf("(S, %d, %d) sub\n", offs[MAT_SUB][s][d], 
                        offs[MAT_SUB][s][d]+diag);
            }

            // open gap (enter DEL)
            if (s-(o+e) >= 0 && d > 0 && 
                           offs[MAT_SUB][s-(o+e)][d-1] != -2 &&
                    diag + offs[MAT_SUB][s-(o+e)][d-1] < truth_len &&
                           offs[MAT_SUB][s-(o+e)][d-1] >= offs[MAT_DEL][s][d]) {
                offs[MAT_DEL][s][d] = offs[MAT_SUB][s-(o+e)][d-1];
                if(print) printf("(D, %d, %d) open\n", offs[MAT_DEL][s][d], 
                        offs[MAT_DEL][s][d]+diag);
            }

            // open gap (enter INS)
            if (s-(o+e) >= 0 && d < mat_len-1 && 
                           offs[MAT_SUB][s-(o+e)][d+1] != -2 &&
                           offs[MAT_SUB][s-(o+e)][d+1]+1 < query_len &&
                    diag + offs[MAT_SUB][s-(o+e)][d+1]+1 < truth_len &&
                    diag + offs[MAT_SUB][s-(o+e)][d+1]+1 >= 0 &&
                           offs[MAT_SUB][s-(o+e)][d+1]+1 >= offs[MAT_INS][s][d]) {
                offs[MAT_INS][s][d] = offs[MAT_SUB][s-(o+e)][d+1]+1;
                if(print) printf("(I, %d, %d) open\n", offs[MAT_INS][s][d], 
                        offs[MAT_INS][s][d]+diag);
            }

            // extend gap (stay DEL)
            if (s-e >= 0 && d > 0 && 
                           offs[MAT_DEL][s-e][d-1] != -2 &&
                    diag + offs[MAT_DEL][s-e][d-1] < truth_len &&
                           offs[MAT_DEL][s-e][d-1] >= offs[MAT_DEL][s][d]) {
                offs[MAT_DEL][s][d] = offs[MAT_DEL][s-e][d-1];
                if(print) printf("(D, %d, %d) extend\n", offs[MAT_DEL][s][d], 
                        offs[MAT_DEL][s][d]+diag);
            }

            // extend gap (stay INS)
            if (s-e >= 0 && d < mat_len-1 && 
                           offs[MAT_INS][s-e][d+1] != -2 &&
                           offs[MAT_INS][s-e][d+1]+1 < query_len &&
                    diag + offs[MAT_INS][s-e][d+1]+1 < truth_len &&
                    diag + offs[MAT_INS][s-e][d+1]+1 >= 0 &&
                           offs[MAT_INS][s-e][d+1]+1 >= offs[MAT_INS][s][d]) {
                offs[MAT_INS][s][d] = offs[MAT_INS][s-e][d+1]+1;
                if(print) printf("(I, %d, %d) extend\n", offs[MAT_INS][s][d], 
                        offs[MAT_INS][s][d]+diag);
            }
        }
    }
}

/******************************************************************************/

void precision_recall_threads_wrapper(
        std::shared_ptr<superclusterData> clusterdata_ptr,
        std::vector< std::vector< std::vector<int> > > sc_groups) {

    if (g.verbosity >= 1 ) 
    for (int i = 0; i < g.thread_nsteps; i++) {
        INFO("  Superclusters using %7.3f to %7.3f GB RAM each (%3d threads): %8d",
                i == 0 ? 0 : g.ram_steps[i-1], g.ram_steps[i], g.thread_steps[i],
                int(sc_groups[i][CTG_IDX].size()));
    }

    int thread_step = g.thread_nsteps-1;
    int start = 0;
    while (thread_step >= 0) {

        // init
        int nthreads = g.thread_steps[thread_step];
        int nscs = sc_groups[thread_step][SC_IDX].size() - start;
        std::vector<std::thread> threads;

        // all threads can solve at least one problem at this level
        if (nscs >= nthreads) {

            // distribute many subproblems evenly among threads
            bool thread2 = thread_step != 0; // max_threads/2+
            for (int t = 0; t < nthreads; t++) {
                int size = nscs / nthreads;
                threads.push_back(std::thread(precision_recall_wrapper,
                            clusterdata_ptr.get(), std::cref(sc_groups),
                            thread_step, start, start+size, thread2, /* print = */ false));
                start += size;
            }
            for (std::thread & t : threads)
                t.join();

            // we fully finished all problems at this size
            if (nscs % nthreads == 0) {
                thread_step--;
                start = 0;
            }

        } else { // not enough work for all threads at this level

            // each thread solves one problem, 
            // at multiple levels until RAM limit reached
            float total_ram = 0;
            while (total_ram < g.max_ram && thread_step >= 0) {
                while (sc_groups[thread_step][SC_IDX].size() == 0) {
                    thread_step--;
                    if (thread_step < 0) break;
                }
                bool thread2 = thread_step != 0; // max_threads/2+
                if (thread_step < 0) break;
                threads.push_back(std::thread(precision_recall_wrapper,
                            clusterdata_ptr.get(), std::cref(sc_groups),
                            thread_step, start, start+1, thread2, /* print = */ false));
                start++;
                total_ram += g.ram_steps[thread_step];
                if (start >= int(sc_groups[thread_step][SC_IDX].size())) {
                    start = 0;
                    thread_step--;
                }
            }

            for (std::thread & t : threads)
                t.join();
        }
    }
}


/******************************************************************************/


void Graph::print() {
    printf("QUERY:\n");
    for (int qn = 0; qn < this->qnodes; qn++) {
        printf("node %d = %s %s\n", qn, this->qseqs[qn].data(), type_strs[this->qtypes[qn]].data());
        printf("\tprev nodes: ");
        for (int n2 : this->qprevs[qn]) printf("%d, ", n2);
        printf("\n");
        printf("\tnext nodes: ");
        for (int n2 : this->qnexts[qn]) printf("%d, ", n2);
        printf("\n");
    }
    printf("\nTRUTH:\n");
    for (int tn = 0; tn < this->tnodes; tn++) {
        printf("node %d = %s %s\n", tn, this->tseqs[tn].data(), type_strs[this->ttypes[tn]].data());
        printf("\tprev nodes: ");
        for (int n2 : this->tprevs[tn]) printf("%d, ", n2);
        printf("\n");
        printf("\tnext nodes: ");
        for (int n2 : this->tnexts[tn]) printf("%d, ", n2);
        printf("\n");
    }
}


/******************************************************************************/


int Graph::get_truth_pos(int truth_node_idx, int truth_idx) {
    int truth_pos = 0;
    for (int tn = 0; tn < truth_node_idx; tn++) {
        truth_pos += this->tseqs[tn].size() - 1; // each tseq starts with '_'
    }
    truth_pos += truth_idx;
    return truth_pos;
}

/******************************************************************************/


/* Build a graph containing all possible paths through query, and selected truth path
   NOTE: start off simple O(V*V) and make sure it's correct (worry about efficiency later)
 */
Graph::Graph(
		std::shared_ptr<ctgSuperclusters> sc, int sc_idx, 
        std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hap) {

    ////////////////////////////////
    // STEP 1: CREATE QUERY NODES //
    ////////////////////////////////

    // initialize helper variables
    std::shared_ptr<ctgVariants> qvars = sc->callset_vars[QUERY];
    std::shared_ptr<ctgVariants> tvars = sc->callset_vars[TRUTH];
    int qvar_beg = std::distance(qvars->superclusters.begin(),
            std::lower_bound(qvars->superclusters.begin(), qvars->superclusters.end(), sc_idx));
    int qvar_end = std::distance(qvars->superclusters.begin(),
            std::upper_bound(qvars->superclusters.begin(), qvars->superclusters.end(), sc_idx));
    int tvar_beg = std::distance(tvars->superclusters.begin(),
            std::lower_bound(tvars->superclusters.begin(), tvars->superclusters.end(), sc_idx));
    int tvar_end = std::distance(tvars->superclusters.begin(),
            std::upper_bound(tvars->superclusters.begin(), tvars->superclusters.end(), sc_idx));
    int ref_beg = sc->get_min_ref_pos(qvar_beg, qvar_end, tvar_beg, tvar_end);
    int ref_end = sc->get_max_ref_pos(qvar_beg, qvar_end, tvar_beg, tvar_end);

    int ref_pos = ref_beg;
    this->sc = sc;
    this->sc_idx = sc_idx;
    // keep track of all node ends. variants are already sorted in order of their starts, but
    // to create the correct reference substrings we need to know where the next end is too
    std::priority_queue<int, std::vector<int>, std::greater<int> > qnode_ends;

    // iterate through all the variants
    this->qnodes = 0;
    for (int var_idx = qvar_beg; var_idx < qvar_end; var_idx++) {
        // skip variants that we've already evaluated
        if (qvars->errtypes[truth_hap][var_idx] != ERRTYPE_UN) { continue; }
        
        // prior to adding the next variant, add as many reference nodes as necessary
        // due to earlier variants that end before it starts
        int var_pos = qvars->poss[var_idx];
        while (ref_pos < var_pos) {
            // get length of this reference segment
            int next_ref_pos; 
            if (!qnode_ends.empty() && qnode_ends.top() < var_pos) { // prev var ends before next starts
                next_ref_pos = qnode_ends.top();
                // multiple vars may end at same location; remove all from node_ends
                while (!qnode_ends.empty() && qnode_ends.top() == next_ref_pos) qnode_ends.pop();
            } else {
                next_ref_pos = var_pos;
            }
            // add reference node
            assert(next_ref_pos <= var_pos);
            this->qnodes++;
            this->qseqs.push_back("_" + ref->fasta.at(ctg).substr(ref_pos, next_ref_pos-ref_pos));
            this->qbegs.push_back(ref_pos - ref_beg);
            this->qends.push_back(next_ref_pos - ref_beg);
            this->qtypes.push_back(TYPE_REF);
            this->qidxs.push_back(-1);
            ref_pos = next_ref_pos;
        }

        // add the variant
        assert(ref_pos == var_pos);
        this->qnodes++;
        this->qseqs.push_back("_" + qvars->alts[var_idx]);
        this->qbegs.push_back(qvars->poss[var_idx] - ref_beg);
        this->qends.push_back(qvars->poss[var_idx] + qvars->rlens[var_idx] - ref_beg);
        this->qtypes.push_back(qvars->types[var_idx]);
        this->qidxs.push_back(var_idx);
        if (qvars->rlens[var_idx]) 
            qnode_ends.push(qvars->poss[var_idx] + qvars->rlens[var_idx]);
    }

    // all variants added, add the remainder of the reference
    while (ref_pos < ref_end) {
        // get the length of this reference segment
        int next_ref_pos; 
        if (!qnode_ends.empty() && qnode_ends.top() < ref_end) { // stop at variant end
            next_ref_pos = qnode_ends.top();
            // multiple vars may end at same location; remove all from node_ends
            while (!qnode_ends.empty() && qnode_ends.top() == next_ref_pos) qnode_ends.pop();
        } else {
            next_ref_pos = ref_end+1; // add remainder of reference in this supercluster
            // TODO: why is +1 necessary for graph compared to generate_ptrs_strs()?
        }
        this->qnodes++;
        this->qseqs.push_back("_" + ref->fasta.at(ctg).substr(ref_pos, next_ref_pos-ref_pos));
        this->qbegs.push_back(ref_pos - ref_beg);
        this->qends.push_back(next_ref_pos - ref_beg);
        this->qtypes.push_back(TYPE_REF);
        this->qidxs.push_back(-1);
        ref_pos = next_ref_pos;
    }

    /////////////////////////////////////
    // STEP 2: SET QUERY NODE POINTERS //
    /////////////////////////////////////
    this->qprevs.resize(this->qnodes);
    this->qnexts.resize(this->qnodes);

    for (int n1 = 0; n1 < this->qnodes; n1++) {
        for (int n2 = 0; n2 < this->qnodes; n2++) {
            if (n2 == n1) continue; // don't compare node to itself
            // add predecessor of current node
            if (this->qbegs[n1] == this->qends[n2] && n2 < n1)
                this->qprevs[n1].push_back(n2);
            // add successor of current node
            if (this->qends[n1] == this->qbegs[n2] && n1 < n2)
                this->qnexts[n1].push_back(n2);
        }
    }

    ////////////////////////////////
    // STEP 3: CREATE TRUTH NODES //
    ////////////////////////////////

    // iterate through all the variants
    ref_pos = ref_beg;
    this->tnodes = 0;
    this->truth = "";
    for (int var_idx = tvar_beg; var_idx < tvar_end; var_idx++) {
        // skip variants that we've already evaluated
        if (tvars->errtypes[truth_hap][var_idx] != ERRTYPE_UN) { continue; }

        if (tvars->var_on_hap(var_idx, truth_hap)) { // ignore variant if on other hap
            // add reference node before next truth variant
            int var_pos = tvars->poss[var_idx];
            if (var_pos > ref_pos) {
                this->tnodes++;
                this->tseqs.push_back("_" + ref->fasta.at(ctg).substr(ref_pos, var_pos-ref_pos));
                this->tbegs.push_back(ref_pos - ref_beg);
                this->ttypes.push_back(TYPE_REF);
                this->tidxs.push_back(-1);
                this->truth += ref->fasta.at(ctg).substr(ref_pos, var_pos-ref_pos);
            }

            // add the truth variant
            this->tnodes++;
            this->tseqs.push_back("_" + tvars->alts[var_idx]);
            this->tbegs.push_back(var_pos - ref_beg);
            this->ttypes.push_back(tvars->types[var_idx]);
            this->tidxs.push_back(var_idx);
            this->truth += tvars->alts[var_idx];
            ref_pos = tvars->poss[var_idx] + tvars->rlens[var_idx];
        }
    }

    // add the remainder of the truth
    this->tnodes++;
    this->tseqs.push_back("_" + ref->fasta.at(ctg).substr(ref_pos, ref_end+1 - ref_pos));
    this->tbegs.push_back(ref_pos - ref_beg);
    this->ttypes.push_back(TYPE_REF);
    this->tidxs.push_back(-1);
    this->truth += ref->fasta.at(ctg).substr(ref_pos, ref_end+1 - ref_pos);

    /////////////////////////////////////
    // STEP 4: SET TRUTH NODE POINTERS //
    /////////////////////////////////////

    this->tprevs.resize(this->tnodes);
    this->tnexts.resize(this->tnodes);
    for (int tni = 0; tni < this->tnodes; tni++) {
        int last_tnode = this->tnodes-1;
        if (tni != last_tnode) {
            this->tnexts[tni].push_back(tni+1);
            this->tprevs[tni+1].push_back(tni);
        }
    }

    this->ref = ref->fasta.at(ctg).substr(ref_beg, ref_end+1 - ref_beg);
}


/******************************************************************************/

void precision_recall_wrapper(
        superclusterData* clusterdata_ptr,
        const std::vector< std::vector< std::vector<int> > > & sc_groups,
        int thread_step, int start, int stop, bool thread2, bool print) {

	// parse sc_idx from grouped superclusters
    if (stop == start) return;
    for (int supclust_idx = start; supclust_idx < stop; supclust_idx++) {
        std::string ctg = clusterdata_ptr->contigs[
            sc_groups[thread_step][CTG_IDX][supclust_idx]];
        int sc_idx = sc_groups[thread_step][SC_IDX][supclust_idx];

        // set superclusters pointer
        std::shared_ptr<ctgSuperclusters> scs = clusterdata_ptr->superclusters[ctg];

        if (print) {
            // print cluster info
            printf("\n\nSupercluster: %d\n", sc_idx);
            for (int c = 0; c < CALLSETS; c++) {
                std::shared_ptr<ctgVariants> vars = scs->callset_vars[c];
                int var_beg = std::distance(vars->superclusters.begin(),
                        std::lower_bound(vars->superclusters.begin(), 
                                         vars->superclusters.end(), sc_idx));
                int var_end = std::distance(vars->superclusters.begin(),
                        std::upper_bound(vars->superclusters.begin(), 
                                         vars->superclusters.end(), sc_idx));
                printf("%s: %d variants (%d-%d)\n", callset_strs[c].data(),
                        var_end-var_beg, var_beg, var_end);
                for (int var_idx = var_beg; var_idx < var_end; var_idx++) {
                    printf("\t\t%s %d\t%s\t%s\t%s\tQ=%f\n", ctg.data(), vars->poss[var_idx], 
                    vars->refs[var_idx].size() ?  vars->refs[var_idx].data() : "_", 
                    vars->alts[var_idx].size() ?  vars->alts[var_idx].data() : "_",
                    gt_strs[vars->orig_gts[var_idx]].data(),
                    vars->var_quals[var_idx]);
                }
            }
        }
        
        // calculate two forward-pass alignments, saving path
        // query1/query2 graph to truth1, query1/query2 graph to truth2
        for (int hi = 0; hi < HAPS; hi++) {
            if (print) printf("HAP %d\n", hi);
            evaluate_variants(scs, sc_idx, clusterdata_ptr->ref, ctg, hi, print);
        }
    }
}


/******************************************************************************/


/* Perform WaveFront Smith-Waterman-Gotoh alignment of two strings, returning 
 * the farthest-reaching reference index of lesser or equal score to that provided.
 * Strings are generated by applying variants to draft ref, and skipping 
 * variants is not allowed. Neither is a diagonal reference transition past 
 * the existing supercluster.
 */
int wf_swg_max_reach(
        const std::string & query, const std::string & truth, 
        std::vector<int> & offs,
        int main_diag, int main_diag_start, int max_score, 
        int x, int o, int e, bool print /* false */, bool reverse /* false */
        ) {

    // init
    int query_len = query.size();
    int truth_len = truth.size();
    int mat_len = query_len + truth_len - 1;
    int main_diag_off = main_diag_start - main_diag;
    int s = 0;
    int s2 = 0;
    int scores = std::max(x, o+e)+1;
    int y = mat_len;
    int z = y * scores;
    offs[MAT_SUB*z + s2*y + query_len-1] = -1;

    while (true) {

        // EXTEND WAVEFRONT (leave INS, DEL forwards)
        if (!reverse) for (int m = MAT_INS; m < MATS; m++) {
            for (int d = 0; d < mat_len; d++) {
                int off = offs[m*z + s2*y + d];
                int diag = d + 1 - query_len;

                if (off >= 0 && off < query_len &&
                        diag+off >= 0 && diag+off < truth_len &&
                        off >= offs[MAT_SUB*z + s2*y + d]) {
                    offs[MAT_SUB*z + s2*y + d] = off;
                    if(print) printf("(S, %d, %d) swap fwd\n", off, off+diag);

                }
            }
        }

        // EXTEND WAVEFRONT (diag, SUB only)
        for (int d = 0; d < mat_len; d++) {
            int off = offs[MAT_SUB*z + s2*y + d];
            int diag = d + 1 - query_len;

            // extend
            while ( (diag != main_diag || off+1 < main_diag_off) &&
                        off != -2 && diag + off >= -1 && 
                        off < query_len - 1 && 
                        diag + off < truth_len - 1) {
                if (query[off+1] == truth[diag+off+1]) off++;
                else break;
            }
            if (off > offs[MAT_SUB*z + s2*y + d])
                if(print) printf("(S, %d, %d) extend\n", off, off+diag);
            offs[MAT_SUB*z + s2*y + d] = off;

            // finish if we've reached the last column
            if (off + diag == truth_len - 1) {
                return truth_len-1;
            }
            if (off == query_len - 1 && off+diag >= 0 && off+diag < truth_len-1)  {
                ERROR("Reached end of query (len %d) before ref (len %d) in max_reach() at (%d, %d)",
                        query_len, truth_len, off, diag+off);
            }

        }
        if (s == max_score) break;

        /* if (print) for (int mi = 0; mi < MATS; mi++) { */
        /*     printf("\n%s matrix\n", type_strs[mi+1].data()); */
        /*     printf("offs %d:", s); */
        /*     for (int di = 0; di < int(query.size() + truth.size()-1); di++) { */
        /*         printf("\t%d", offs[mi*z + s2*y + di]); */
        /*     } */
        /*     printf("\n"); */
        /* } */

        // NEXT WAVEFRONT
        s++; s2++;
        if (s2 == scores) s2 = 0;
        // init new row to all invalid
        for (int m = MAT_INS; m < MATS; m++) {
            for (int d = 0; d < mat_len; d++) {
                offs[m*z + s2*y + d] = -2;
            }
        }
        if (print) printf("\nscore = %d\n", s);

        for (int d = 0; d < mat_len; d++) {
            int diag = d + 1 - query_len;

            // sub (in SUB)
            int p = s - x;
            int p2 = s2 - x;
            if (p2 < 0) p2 += scores;
            if (p >= 0 && offs[MAT_SUB*z + p2*y + d] != -2 && 
                            offs[MAT_SUB*z + p2*y + d]+1 < query_len &&
                     diag + offs[MAT_SUB*z + p2*y + d]+1 < truth_len &&
                            offs[MAT_SUB*z + p2*y + d]+1 >= offs[MAT_SUB*z + s2*y + d]) {
                offs[MAT_SUB*z + s2*y + d] = offs[MAT_SUB*z + p2*y + d] + 1;
                if(print) printf("(S, %d, %d) sub\n", offs[MAT_SUB*z + s2*y + d], 
                        offs[MAT_SUB*z + s2*y + d]+diag);
            }

            // open gap (enter DEL, open fwd only)
            p = reverse ? s - e : s - (o+e);
            p2 = reverse ? s2 - e : s2 - (o+e);
            if (p2 < 0) p2 += scores;
            if (p >= 0 && d > 0 && 
                           offs[MAT_SUB*z + p2*y + d-1] != -2 &&
                    diag + offs[MAT_SUB*z + p2*y + d-1] < truth_len &&
                           offs[MAT_SUB*z + p2*y + d-1] >= offs[MAT_DEL*z + s2*y + d]) {
                offs[MAT_DEL*z + s2*y + d] = offs[MAT_SUB*z + p2*y + d-1];
                if(print) printf("(D, %d, %d) open\n", offs[MAT_DEL*z + s2*y + d], 
                        offs[MAT_DEL*z + s2*y + d]+diag);
            }
            // open gap (enter INS, open fwd only)
            if (p >= 0 && d < mat_len-1 && 
                           offs[MAT_SUB*z + p2*y + d+1] != -2 &&
                           offs[MAT_SUB*z + p2*y + d+1]+1 < query_len &&
                    diag + offs[MAT_SUB*z + p2*y + d+1]+1 < truth_len &&
                    diag + offs[MAT_SUB*z + p2*y + d+1]+1 >= 0 &&
                           offs[MAT_SUB*z + p2*y + d+1]+1 >= offs[MAT_INS*z + s2*y + d]) {
                offs[MAT_INS*z + s2*y + d] = offs[MAT_SUB*z + p2*y + d+1]+1;
                if(print) printf("(I, %d, %d) open\n", offs[MAT_INS*z + s2*y + d], 
                        offs[MAT_INS*z + s2*y + d]+diag);
            }

            // leave INDEL (open rev only)
            p = s - o;
            p2 = s2 - o;
            if (p2 < 0) p2 += scores;
            if (reverse && p >= 0) {
                for (int m = MAT_INS; m < MATS; m++) {
                    if (        offs[m*z + p2*y + d] >= 0 && 
                                offs[m*z + p2*y + d] < query_len &&
                         diag + offs[m*z + p2*y + d] >= 0 && 
                         diag + offs[m*z + p2*y + d] < truth_len &&
                                offs[m*z + p2*y + d] > offs[MAT_SUB*z + s2*y + d]) {
                        offs[MAT_SUB*z + s2*y + d] = offs[m*z + p2*y + d];
                        if(print) printf("(S, %d, %d) swap rev\n", 
                                offs[m*z + p2*y + d], diag+offs[m*z + p2*y + d]);
                    }
                }
            }

            // extend gap (stay DEL)
            p = s - e;
            p2 = s2 - e;
            if (p2 < 0) p2 += scores;
            if (p >= 0 && d > 0 && 
                           offs[MAT_DEL*z + p2*y + d-1] != -2 &&
                    diag + offs[MAT_DEL*z + p2*y + d-1] < truth_len &&
                           offs[MAT_DEL*z + p2*y + d-1] >= offs[MAT_DEL*z + s2*y + d]) {
                offs[MAT_DEL*z + s2*y + d] = offs[MAT_DEL*z + p2*y + d-1];
                if(print) printf("(D, %d, %d) extend\n", offs[MAT_DEL*z + s2*y + d], 
                        offs[MAT_DEL*z + s2*y + d]+diag);
            }
            // extend gap (stay INS)
            if (p >= 0 && d < mat_len-1 && 
                           offs[MAT_INS*z + p2*y + d+1] != -2 &&
                           offs[MAT_INS*z + p2*y + d+1]+1 < query_len &&
                    diag + offs[MAT_INS*z + p2*y + d+1]+1 < truth_len &&
                    diag + offs[MAT_INS*z + p2*y + d+1]+1 >= 0 &&
                           offs[MAT_INS*z + p2*y + d+1]+1 >= offs[MAT_INS*z + s2*y + d]) {
                offs[MAT_INS*z + s2*y + d] = offs[MAT_INS*z + p2*y + d+1]+1;
                if(print) printf("(I, %d, %d) extend\n", offs[MAT_INS*z + s2*y + d], 
                        offs[MAT_INS*z + s2*y + d]+diag);
            }
        }
    } // end reach

    // get max reach
    int max_reach = 0;
    for (int s2 = 0; s2 < scores; s2++) {
        for (int m = 0; m < MATS; m++) {
            for (int d = 0; d < mat_len; d++) {
                int off = offs[m*z + s2*y + d];
                int diag = d + 1 - query_len;
                if (off >= 0 && off < query_len &&
                        diag+off >= 0 && diag+off < truth_len)
                    max_reach = std::max(max_reach, diag + off);
            }
        }
    }
    return max_reach;
}


/******************************************************************************/


void wf_ed(const std::string & query, const std::string & truth, int & s, bool print) {

    // alignment
    int query_len = query.size();
    int truth_len = truth.size();

    // early exit if either string is empty
    if (!query_len) { s = truth_len; return; }
    if (!truth_len) { s = query_len; return; }
    s = 0;

    int mat_len = query_len + truth_len - 1;
    std::vector< std::vector<int> > offs;
    offs.push_back(std::vector<int>(mat_len, -2));
    offs[0][query_len-1] = -1;
    bool done = false;
    while (true) {

        // EXTEND WAVEFRONT
        for (int d = 0; d < mat_len; d++) {
            int off = offs[s][d];
            int diag = d + 1 - query_len;

            // don't allow starting from untouched cells
            if (off == -2) continue;

            // check that it's within matrix
            if (diag + off + 1 < 0) continue;
            if (off > query_len - 1) continue;
            if (diag + off > truth_len - 1) continue;

            // extend
            while (off < query_len - 1 &&
                   diag + off < truth_len - 1) {
                if (query[off+1] == truth[diag+off+1]) off++;
                else break;
            }
            if (off > offs[s][d])
                if(print) printf("(%d, %d) extend\n", off, off+diag);
            offs[s][d] = off;

            // finish if done
            if (off == query_len - 1 && off + diag == truth_len - 1)
            { done = true; break; }

        }
        if (done) break;

        // debug print
        if(print) printf("offs %d:", s);
        for (int di = 0; di < int(query.size() + truth.size()-1); di++) {
            if(print) printf("\t%d", offs[s][di]);
        }
        if(print) printf("\n");

        // NEXT WAVEFRONT
        offs.push_back(std::vector<int>(mat_len, -2));
        s++;
        if(print) printf("\nscore = %d\n", s);
        for (int d = 0; d < mat_len; d++) {
            int diag = d + 1 - query_len;

            // SUB
            if (s-1 >= 0 && offs[s-1][d] != -2 &&
                            offs[s-1][d]+1 < query_len &&
                     diag + offs[s-1][d]+1 < truth_len &&
                            offs[s-1][d]+1 >= offs[s][d]) {
                offs[s][d] = offs[s-1][d] + 1;
                if(print) printf("(%d, %d) sub\n", offs[s][d], offs[s][d]+diag);
            }

            // DEL
            if (s-1 >= 0 && d > 0 &&
                           offs[s-1][d-1] != -2 &&
                    diag + offs[s-1][d-1] < truth_len &&
                           offs[s-1][d-1] >= offs[s][d]) {
                offs[s][d] = offs[s-1][d-1];
                if(print) printf("(%d, %d) del\n", offs[s][d], offs[s][d]+diag);
            }

            // INS
            if (s-1 >= 0 && d < mat_len-1 &&
                           offs[s-1][d+1] != -2 &&
                           offs[s-1][d+1]+1 < query_len &&
                    diag + offs[s-1][d+1]+1 < truth_len &&
                    diag + offs[s-1][d+1]+1 >= -1 &&
                           offs[s-1][d+1]+1 >= offs[s][d]) {
                offs[s][d] = offs[s-1][d+1]+1;
                if(print) printf("(%d, %d) ins\n", offs[s][d], offs[s][d]+diag);
            }
        }
    }
}

