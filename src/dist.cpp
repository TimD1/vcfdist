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


/* Extract the FP and FN variants from evaluated superclusters into new variant lists. 
 * These variants will be re-evaluated in combination with larger variants that skipped
 * the first round of evaluation.
 */
void extract_errors(std::shared_ptr<superclusterData> sc_data_ptr,
        std::shared_ptr<variantData> query_fp_ptr,
        std::shared_ptr<variantData> truth_fn_ptr) {

    // initialize query variant data
    query_fp_ptr->ref = sc_data_ptr->ref;
    query_fp_ptr->callset = QUERY;
    query_fp_ptr->contigs = sc_data_ptr->contigs;
    query_fp_ptr->lengths = sc_data_ptr->lengths;
    query_fp_ptr->ploidy = sc_data_ptr->ploidy;
    query_fp_ptr->sample = sc_data_ptr->samples[QUERY];
    query_fp_ptr->filename = sc_data_ptr->filenames[QUERY];

    // initialize truth variant data
    truth_fn_ptr->ref = sc_data_ptr->ref;
    truth_fn_ptr->callset = TRUTH;
    truth_fn_ptr->contigs = sc_data_ptr->contigs;
    truth_fn_ptr->lengths = sc_data_ptr->lengths;
    truth_fn_ptr->ploidy = sc_data_ptr->ploidy;
    truth_fn_ptr->sample = sc_data_ptr->samples[TRUTH];
    truth_fn_ptr->filename = sc_data_ptr->filenames[TRUTH];

    for (const std::string & ctg : sc_data_ptr->contigs) {
        // create empty query contig variants
        for (int hi = 0; hi < HAPS; hi++) {
            std::shared_ptr<ctgVariants> ctg_vars(new ctgVariants(ctg));
            query_fp_ptr->variants[hi][ctg] = ctg_vars;
        }

        // copy FP query variants over to new list
        std::shared_ptr<ctgVariants> qvars = sc_data_ptr->superclusters[ctg]->callset_vars[QUERY];
        std::vector<int> qvars_to_remove;
        for (int qi = 0; qi < qvars->n; qi++) {
            // TODO: for now, only pull out variants that are FP on both haps
            if (qvars->errtypes[HAP1][qi] == ERRTYPE_FP && qvars->errtypes[HAP2][qi] == ERRTYPE_FP) {
                qvars_to_remove.push_back(qi);
                for (int hi = 0; hi < HAPS; hi++) {
                    if (qvars->var_on_hap(qi, hi)) {
                        query_fp_ptr->variants[hi][ctg]->add_var(qvars, qi);
                    }
                }
            }
        }

        // remove FP query variants from evaluated variants
        qvars->remove_vars(qvars_to_remove);
        
        // create empty truth contig variants
        for (int hi = 0; hi < HAPS; hi++) {
            std::shared_ptr<ctgVariants> ctg_vars(new ctgVariants(ctg));
            truth_fn_ptr->variants[hi][ctg] = ctg_vars;
        }

        // copy FN truth variants over to new list
        std::shared_ptr<ctgVariants> tvars = sc_data_ptr->superclusters[ctg]->callset_vars[TRUTH];
        std::vector<int> tvars_to_remove;
        for (int ti = 0; ti < tvars->n; ti++) {
            // TODO: for now, only pull out variants that are FN on both haps
            if (tvars->errtypes[HAP1][ti] == ERRTYPE_FN && tvars->errtypes[HAP2][ti] == ERRTYPE_FN) {
                tvars_to_remove.push_back(ti);
                for (int hi = 0; hi < HAPS; hi++) {
                    if (tvars->var_on_hap(ti, hi)) {
                        truth_fn_ptr->variants[hi][ctg]->add_var(tvars, ti);
                    }
                }
            }
        }

        // remove FN truth variants from evaluated variants
        tvars->remove_vars(tvars_to_remove);
    }
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


/* This function is an exact copy of wfa_calc_prec_recall_aln(), but without ptrs, and the 'score'
   dimension of offs is size 2. This function should be used when backtracking is not necessary, and
   uses less memory.
 */
int gwfa_aln(
        const std::shared_ptr<Graph> graph,
        std::vector< std::vector< std::vector< std::vector<int> > > > & offs,
        bool print
        ) {

    // init: offs[query_node][truth_node][score][diagonal] = offset
    int s = 0;
    bool done = false;
    for (int qni = 0; qni < graph->qnodes; qni++) {
        offs.push_back(std::vector< std::vector< std::vector<int> > >(graph->tnodes));
        for (int tni = 0; tni < graph->tnodes; tni++) {
            int mat_len = graph->qseqs[qni].length() + graph->tseqs[tni].length() - 1;
            offs[qni][tni].push_back(std::vector<int>(mat_len, -2));
            offs[qni][tni].push_back(std::vector<int>(mat_len, -2));
        }
    }
    offs[0][0][s][graph->qseqs[0].length()-1] = -1; // main diag

    while (s <= g.max_dist) {

        // EXTEND WAVEFRONT
        // NOTE: because the graph is a topologically sorted DAG, we can safely assume that all
        // extensions that reach into other matrices will always have a larger qni/tni, and be 
        // fully extended later in this double for loop
        for (int qni = 0; qni < graph->qnodes; qni++) {
            for (int tni = 0; tni < graph->tnodes; tni++) {
                int query_len = graph->qseqs[qni].length();
                int truth_len = graph->tseqs[tni].length();
                int mat_len = query_len + truth_len - 1;
                for (int d = 0; d < mat_len; d++) {
                    int off = offs[qni][tni][s&1][d];
                    int diag = d + 1 - query_len;

                    // extend
                    while (off != -2 && diag + off >= -1 && 
                           off < query_len - 1 && 
                           diag + off < truth_len - 1) {
                        if (graph->qseqs[qni][off+1] == graph->tseqs[tni][diag+off+1]) off++;
                        else break;
                    }
                    if (off == query_len - 1) { // extend into next query matrices
                        for (int next_qni : graph->qnexts[qni]) {
                            int next_query_len = graph->qseqs[next_qni].length();
                            int next_d = d + next_query_len - 1;  // derived
                            if (offs[next_qni][tni][s&1][next_d] == -2) {
                                offs[next_qni][tni][s&1][next_d] = 0;
                            }
                        }
                    }
                    if (diag + off == truth_len - 1) {
                        for (int next_tni : graph->tnexts[tni]) {
                            int next_query_len = graph->qseqs[qni].length();
                            int next_d = next_query_len - off - 1;
                            if (offs[qni][next_tni][s&1][next_d] == -2) {
                                offs[qni][next_tni][s&1][next_d] = off;
                            }
                        }
                    }

                    // debug printing
                    if (print) {
                        if (off > offs[qni][tni][s&1][d]) {
                            printf("(%d, %d, %d, %d) extend\n", qni, tni, off, off+diag);
                        }
                        if (off == query_len - 1) {
                            for (int next_qni : graph->qnexts[qni]) {
                                printf("(%d, %d, 0, %d) query extend\n", next_qni, tni, d);
                            }
                        }
                        if (diag + off == truth_len - 1) {
                            for (int next_tni : graph->tnexts[tni]) {
                                printf("(%d, %d, %d, 0) truth extend\n", qni, next_tni, off);
                            }
                        }
                    }

                    // update offset
                    offs[qni][tni][s&1][d] = off;

                    // finish if done
                    if (qni == graph->qnodes - 1 && 
                        tni == graph->tnodes - 1 &&
                        off == query_len - 1 && 
                        off + diag == truth_len - 1)
                    { done = true; break; }
                }
            }
        }
        if (done) break;

        // debug print
        if (print)
        for (int qni = 0; qni < graph->qnodes; qni++) {
            for (int tni = 0; tni < graph->tnodes; tni++) {
                printf("\n(%d, %d) matrix\n", qni, tni);
                printf("offs %d:", s);
                int query_len = graph->qseqs[qni].length();
                int truth_len = graph->tseqs[tni].length();
                int mat_len = query_len + truth_len - 1;
                for (int di = 0; di < mat_len; di++) {
                    printf("\t%d", offs[qni][tni][s&1][di]);
                }
                printf("\n");
            }
        }

        // NEXT WAVEFRONT
        s++;
        if(print) printf("\nscore = %d\n", s);

        for (int qni = 0; qni < graph->qnodes; qni++) {
            for (int tni = 0; tni < graph->tnodes; tni++) {
                int query_len = graph->qseqs[qni].length();
                int truth_len = graph->tseqs[tni].length();
                int mat_len = query_len + truth_len - 1;
                for (int d = 0; d < mat_len; d++) {
                    int diag = d + 1 - query_len;

                    // substitution
                    if (       offs[qni][tni][(s-1)&1][d] != -2 && 
                               offs[qni][tni][(s-1)&1][d]+1 < query_len &&
                        diag + offs[qni][tni][(s-1)&1][d]+1 < truth_len &&
                               offs[qni][tni][(s-1)&1][d]+1 > offs[qni][tni][s&1][d]) {
                        offs[qni][tni][s&1][d] = offs[qni][tni][(s-1)&1][d] + 1;
                        if(print) printf("(%d, %d, %d, %d) sub\n", qni, tni, offs[qni][tni][s&1][d], 
                                offs[qni][tni][s&1][d]+diag);
                    }

                    // deletion
                    if (d > 0 && 
                               offs[qni][tni][(s-1)&1][d-1] != -2 &&
                        diag + offs[qni][tni][(s-1)&1][d-1] < truth_len &&
                               offs[qni][tni][(s-1)&1][d-1] > offs[qni][tni][s&1][d]) {
                        offs[qni][tni][s&1][d] = offs[qni][tni][(s-1)&1][d-1];
                        if(print) printf("(%d, %d, %d, %d) del\n", qni, tni, offs[qni][tni][s&1][d], 
                                offs[qni][tni][s&1][d]+diag);
                    }

                    // insertion
                    if (d < mat_len-1 && 
                               offs[qni][tni][(s-1)&1][d+1] != -2 &&
                               offs[qni][tni][(s-1)&1][d+1]+1 < query_len &&
                        diag + offs[qni][tni][(s-1)&1][d+1]+1 < truth_len &&
                        diag + offs[qni][tni][(s-1)&1][d+1]+1 >= 0 &&
                               offs[qni][tni][(s-1)&1][d+1]+1 > offs[qni][tni][s&1][d]) {
                        offs[qni][tni][s&1][d] = offs[qni][tni][(s-1)&1][d+1]+1;
                        if(print) printf("(%d, %d, %d, %d) ins\n", qni, tni, offs[qni][tni][s&1][d], 
                                offs[qni][tni][s&1][d]+diag);
                    }
                }
            }
        }
    }
    return s;
}


/******************************************************************************/


int wfa_calc_prec_recall_aln(
        const std::shared_ptr<Graph> graph,
        std::vector<uint32_t> & ptrs,
        std::vector<int> & offs,
        std::unordered_map<idx3, uint64_t> & node_offs,
        bool print
        ) {

    // fail if the graph contains too many nodes
    if (graph->tnodes >= int(1 << NODE_BITS)) {
        ERROR("Alignment graph for supercluster exceeds maximum truth nodes (%d > %d)",
                graph->tnodes, int(1 << NODE_BITS)-1);
    }
    if (graph->qnodes >= int(1 << NODE_BITS)) {
        ERROR("Alignment graph for supercluster exceeds maximum query nodes (%d > %d)",
                graph->qnodes, int(1 << NODE_BITS)-1);
    }

    // init: offs/ptrs [ node_offs[(qnode, tnode, score)] + idx ] = offset/pointer
    int s = 0;
    bool done = false;
    int init_query_len = graph->qseqs[0].length();
    int init_truth_len = graph->tseqs[0].length();
    int init_mat_len = init_query_len + init_truth_len - 1;
    node_offs[idx3(0, 0, 0)] = 0;
    ptrs.insert(ptrs.end(), init_mat_len, PTR_NONE);
    offs.insert(offs.end(), init_mat_len, -2);
    offs[init_query_len-1] = -1; // main diag
    ptrs[init_query_len-1] = PTR_MAT;
    if (print) printf("{%d, %d} score %d = %lu-%lu\n", 0, 0, s,
            node_offs[idx3(0, 0, 0)], ptrs.size());

    while (s <= g.max_dist) {

        // EXTEND WAVEFRONT
        // NOTE: because the graph is a topologically sorted DAG, we can safely assume that all
        // extensions that reach into other matrices will always have a larger qni/tni, and be 
        // fully extended later in this double for loop
        for (int qni = 0; qni < graph->qnodes; qni++) {
            for (int tni = 0; tni < graph->tnodes; tni++) {
                if (print) printf("processing node {%d, %d}\n", qni, tni);
                idx3 node(qni, tni, s);
                if (!contains(node_offs, node)) {
                    continue;
                }
                int query_len = graph->qseqs[qni].length();
                int truth_len = graph->tseqs[tni].length();
                int mat_len = query_len + truth_len - 1;
                for (int d = 0; d < mat_len; d++) {
                    int off = offs[node_offs[node] + d];
                    int diag = d + 1 - query_len;

                    // extend within matrix
                    while (off != -2 && diag + off >= -1 && 
                           off < query_len - 1 && 
                           diag + off < truth_len - 1) {
                        if (graph->qseqs[qni][off+1] == graph->tseqs[tni][diag+off+1]) {
                            off++;
                            if (print) printf("(%d, %d, %d, %d) extend\n", qni, tni, off, off+diag);
                        }
                        else break;
                    }

                    // extend into next query and truth matrices
                    if (off == query_len - 1 and diag + off == truth_len - 1) {
                        for (int next_qni : graph->qnexts[qni]) {
                            for (int next_tni : graph->tnexts[tni]) {
                                if (print) printf("(%d, %d, 0, 0) query and truth extend\n", 
                                            next_qni, next_tni);

                                // allocate memory for new node
                                idx3 next_node(next_qni, next_tni, s);
                                int next_query_len = graph->qseqs[next_qni].length();
                                int next_truth_len = graph->tseqs[next_tni].length();
                                int next_mat_len = next_query_len + next_truth_len - 1;
                                if (!contains(node_offs, next_node)) {
                                    node_offs[next_node] = ptrs.size();
                                    ptrs.insert(ptrs.end(), next_mat_len, PTR_NONE);
                                    offs.insert(offs.end(), next_mat_len, -2);
                                    if (print) printf("{%d, %d} score %d = %lu-%lu\n", next_qni,
                                            next_tni, s, node_offs[next_node], ptrs.size());
                                }
                                // align into new node
                                int next_d = next_query_len - 1;
                                if (offs[node_offs[next_node] + next_d] == -2) {
                                    offs[node_offs[next_node] + next_d] = 0;
                                    ptrs[node_offs[next_node] + next_d] = PTR_MAT |
                                        (qni << (PTR_BITS + NODE_BITS)) | (tni << PTR_BITS);
                                }
                            }
                        }
                    }

                    // extend into next query matrices
                    if (off == query_len - 1) {
                        for (int next_qni : graph->qnexts[qni]) {
                            if (print) printf("(%d, %d, 0, %d) query extend\n", next_qni, tni, d);

                            // allocate memory for new node
                            idx3 next_node(next_qni, tni, s);
                            int next_query_len = graph->qseqs[next_qni].length();
                            int next_mat_len = next_query_len + truth_len - 1;
                            if (!contains(node_offs, next_node)) {
                                node_offs[next_node] = ptrs.size();
                                ptrs.insert(ptrs.end(), next_mat_len, PTR_NONE);
                                offs.insert(offs.end(), next_mat_len, -2);
                                if (print) printf("{%d, %d} score %d = %lu-%lu\n", next_qni, 
                                        tni, s, node_offs[next_node], ptrs.size());
                            }

                            // align into new node
                            // TODO: segfault occurs here!
                            int next_d = d + next_query_len - 1;  // derived
                            if (offs[node_offs[next_node] + next_d] == -2) {
                                offs[node_offs[next_node] + next_d] = 0;
                                // PTR_MAT | PTR_INS means vertical match
                                ptrs[node_offs[next_node] + next_d] = PTR_MAT | PTR_INS |
                                    (qni << (PTR_BITS + NODE_BITS)) | (tni << PTR_BITS);
                            }
                        }
                    }

                    // extend into next truth matrices
                    if (diag + off == truth_len - 1) {
                        for (int next_tni : graph->tnexts[tni]) {
                            if (print) printf("(%d, %d, %d, 0) truth extend\n", qni, next_tni, off);

                            // allocate memory for new node
                            idx3 next_node(qni, next_tni, s);
                            int next_truth_len = graph->tseqs[next_tni].length();
                            int next_mat_len = query_len + next_truth_len - 1;
                            if (!contains(node_offs, next_node)) {
                                node_offs[next_node] = ptrs.size();
                                ptrs.insert(ptrs.end(), next_mat_len, PTR_NONE);
                                offs.insert(offs.end(), next_mat_len, -2);
                                if (print) printf("{%d, %d} score %d = %lu-%lu\n", qni, next_tni, s,
                                        node_offs[next_node], ptrs.size());
                            }

                            // align into new node
                            int next_d = query_len - off - 1;
                            if (offs[node_offs[next_node] + next_d] == -2) {
                                offs[node_offs[next_node] + next_d] = off;
                                // PTR_MAT | PTR_DEL means horizontal match
                                ptrs[node_offs[next_node] + next_d] = PTR_MAT | PTR_DEL |
                                    (qni << (PTR_BITS + NODE_BITS)) | (tni << PTR_BITS);
                            }
                        }
                    }

                    // update offset
                    offs[node_offs[node] + d] = off;

                    // finish if done
                    if (qni == graph->qnodes - 1 && 
                        tni == graph->tnodes - 1 &&
                        off == query_len - 1 && 
                        off + diag == truth_len - 1)
                    { done = true; break; }
                }
            }
        }
        if (done) break;

        // debug print
        if (print)
        for (int qni = 0; qni < graph->qnodes; qni++) {
            for (int tni = 0; tni < graph->tnodes; tni++) {
                idx3 node(qni, tni, s);
                printf("\n(%d, %d) matrix\n", qni, tni);
                printf("contains(%d, %d, %d) = %s\n", qni, tni, s, contains(node_offs, node) ? "true" : "false");
                if (contains(node_offs, node)) {
                    printf("node offset: %lu\n", node_offs[node]);
                    printf("offs %d:", s);
                    int query_len = graph->qseqs[qni].length();
                    int truth_len = graph->tseqs[tni].length();
                    int mat_len = query_len + truth_len - 1;
                    for (int di = 0; di < mat_len; di++) {
                        printf("\t%d", offs[node_offs[node] + di]);
                    }
                    printf("\n");
                } else {
                    printf("unallocated\n");
                }
            }
        }

        // NEXT WAVEFRONT
        s++;
        for (int qni = 0; qni < graph->qnodes; qni++) {
            for (int tni = 0; tni < graph->tnodes; tni++) {

                // for each prev node we're aligning in, add new node for updated score
                idx3 prev_node(qni, tni, s-1);
                idx3 node(qni, tni, s);
                if (!contains(node_offs, prev_node)) {
                    continue;
                }
                int query_len = graph->qseqs[qni].length();
                int truth_len = graph->tseqs[tni].length();
                int mat_len = query_len + truth_len - 1;
                node_offs[node] = ptrs.size();
                ptrs.insert(ptrs.end(), mat_len, PTR_NONE);
                offs.insert(offs.end(), mat_len, -2);

                for (int d = 0; d < mat_len; d++) {
                    int diag = d + 1 - query_len;

                    // substitution
                    if (       offs[node_offs[prev_node] + d] != -2 && 
                               offs[node_offs[prev_node] + d]+1 < query_len &&
                        diag + offs[node_offs[prev_node] + d]+1 < truth_len &&
                               offs[node_offs[prev_node] + d]+1 > offs[node_offs[node] + d]) {
                        offs[node_offs[node] + d] = offs[node_offs[prev_node] + d] + 1;
                        ptrs[node_offs[node] + d] = PTR_SUB;
                        if(print) printf("(%d, %d, %d, %d) sub\n", qni, tni, offs[node_offs[node] + d], 
                                offs[node_offs[node] + d]+diag);
                    }

                    // deletion
                    if (d > 0 && 
                               offs[node_offs[prev_node] + d-1] != -2 &&
                        diag + offs[node_offs[prev_node] + d-1] < truth_len &&
                               offs[node_offs[prev_node] + d-1] > offs[node_offs[node] + d]) {
                        offs[node_offs[node] + d] = offs[node_offs[prev_node] + d-1];
                        ptrs[node_offs[node] + d] = PTR_DEL;
                        if(print) printf("(%d, %d, %d, %d) del\n", qni, tni, offs[node_offs[node] + d], 
                                offs[node_offs[node] + d]+diag);
                    }

                    // insertion
                    if (d < mat_len-1 && 
                               offs[node_offs[prev_node] + d+1] != -2 &&
                               offs[node_offs[prev_node] + d+1]+1 < query_len &&
                        diag + offs[node_offs[prev_node] + d+1]+1 < truth_len &&
                        diag + offs[node_offs[prev_node] + d+1]+1 >= 0 &&
                               offs[node_offs[prev_node] + d+1]+1 > offs[node_offs[node] + d]) {
                        offs[node_offs[node] + d] = offs[node_offs[prev_node] + d+1]+1;
                        ptrs[node_offs[node] + d] = PTR_INS;
                        if(print) printf("(%d, %d, %d, %d) ins\n", qni, tni, offs[node_offs[node] + d], 
                                offs[node_offs[node] + d]+diag);
                    }
                }
            }
        }
    }
    return s;
}

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
            idx4 y(x.qni, x.tni, x.qi+1, x.ti+1);
            if (y.qi < int(graph->qseqs[y.qni].length()) && y.ti < int(graph->tseqs[y.tni].length()) &&
                    graph->qseqs[y.qni][y.qi] == graph->tseqs[y.tni][y.ti]) {
                if (!contains(done, y) && !contains(curr_wave, y)) {
                    /* if (print) printf("      y = node (%d, %d) cell (%d, %d)\n", y.qni, y.tni, y.qi, y.ti); */
                    queue.push(y); curr_wave.insert(y); ptrs[y] = x;
                }
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


std::vector<idx4> parse_wfa_path(const std::shared_ptr<Graph> graph, int s,
        const std::vector<uint32_t> & ptrs,
        const std::vector<int> & offs,
        const std::unordered_map<idx3, uint64_t> & node_offs, bool print) {

    std::vector<idx4> path;
    int qni = graph->qnodes-1;
    int tni = graph->tnodes-1;
    int diag = graph->tseqs[tni].length() - graph->qseqs[qni].length();
    int d = diag + graph->qseqs[qni].length() - 1;
    idx3 init_node(qni, tni, s);
    int off = offs[node_offs.at(init_node) + d];
    uint32_t ptr = ptrs[node_offs.at(init_node) + d];
    path.push_back(idx4(qni, tni, off, diag+off));
    if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);

    while (qni > 0 || tni > 0 || diag+off > 0) {
        if (print) printf("at {%d, %d} (%d, %d) ptr [%u|%u|%u]\n",
                qni, tni, off, diag + off, (ptr & QNODE_MASK) >> (PTR_BITS + NODE_BITS),
                (ptr & TNODE_MASK) >> PTR_BITS, ptr & PTR_MASK);

        // pointer movement into other matrix
        if (ptr & PTR_MAT) {
            // move into new query matrix (hit top row)
            if (ptr & PTR_INS) {
                if (print) printf("insertion, new query matrix\n");
                assert(diag >= 0);
                while (off > 0) {
                    off--;
                    path.push_back(idx4(qni, tni, off, diag+off));
                    if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);
                }
                int prev_qni = (ptr & QNODE_MASK) >> (PTR_BITS + NODE_BITS);
                if (print) printf("hit top row\n");
                path.push_back(idx4(prev_qni, tni, graph->qseqs[prev_qni].length()-1, diag));
                if (print) printf("path {%d, %d} (%d, %d)\n", prev_qni, tni, 
                        int(graph->qseqs[prev_qni].length())-1, diag);
                qni = prev_qni;
                diag = diag - graph->qseqs[qni].length() + 1;
                d = diag + graph->qseqs[qni].length()-1;

            // move into new truth matrix (hit left col)
            } else if (ptr & PTR_DEL) {
                if (print) printf("deletion, new truth matrix\n");
                assert(diag <= 0);
                while (diag + off > 0) {
                    off--;
                    path.push_back(idx4(qni, tni, off, diag+off));
                        if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);
                }
                int prev_tni = (ptr & TNODE_MASK) >> PTR_BITS;
                if (print) printf("hit left col\n");
                path.push_back(idx4(qni, prev_tni, -diag, graph->tseqs[prev_tni].length()-1));
                if (print) printf("path {%d, %d} (%d, %d)\n", qni, prev_tni, -diag,
                        int(graph->tseqs[prev_tni].length())-1);
                tni = prev_tni;
                diag = diag + graph->tseqs[prev_tni].length()-1;
                d = diag + graph->qseqs[qni].length()-1;
            }

            // check if final backtracking matrix, main diagonal
            else if (qni == 0 && tni == 0 && diag == 0) { 
                if (print) printf("final matrix\n");
                while (--off >= -1) {
                    path.push_back(idx4(qni, tni, off, diag+off));
                    if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);
                }
                break;
            }

            // move into new truth and query matrix (main diagonal, exact top left corner)
            else if (diag == 0) {
                if (print) printf("match, new truth and query matrices\n");
                while (off > 0) {
                    off--;
                    path.push_back(idx4(qni, tni, off, diag+off));
                    if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);
                }
                int prev_qni = (ptr & QNODE_MASK) >> (PTR_BITS + NODE_BITS);
                int prev_tni = (ptr & TNODE_MASK) >> PTR_BITS;
                if (print) printf("hit top left corner\n");
                path.push_back(idx4(prev_qni, prev_tni, 
                            graph->qseqs[prev_qni].length()-1,
                            graph->tseqs[prev_tni].length()-1));
                if (print) printf("path {%d, %d} (%d, %d)\n", prev_qni, prev_tni, 
                        int(graph->qseqs[prev_qni].length())-1,
                        int(graph->tseqs[prev_tni].length())-1);
                qni = prev_qni;
                tni = prev_tni;
                diag = graph->tseqs[tni].length() - graph->qseqs[qni].length();
                d = diag + graph->qseqs[qni].length()-1;
            }

            else {
                ERROR("Unexpected pointer value at {%d, %d} (%d, %d) "
                        "(prev=(%u, %u), ptr=%u) during WFA backtrack.", 
                        qni, tni, off, diag + off,
                        (ptr & QNODE_MASK) >> (PTR_BITS + NODE_BITS),
                        (ptr & TNODE_MASK) >> PTR_BITS,
                        ptr & PTR_MASK);
            }
        }

        else {  // same matrix
            idx3 node(qni, tni, s);
            idx3 prev_node(qni, tni, s-1);
            if (print) printf("same matrix\n");
            int prev_off = 0;
            switch (ptrs[node_offs.at(node) + d] & PTR_MASK) {
                case PTR_SUB:
                    if (print) printf("substitution\n");
                    prev_off = offs[node_offs.at(prev_node) + d];
                    s--;
                    while (--off >= prev_off) {
                        path.push_back(idx4(qni, tni, off, diag+off));
                        if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);
                    }
                    break;
                case PTR_INS:
                    if (print) printf("insertion\n");
                    prev_off = offs[node_offs.at(prev_node) + d+1];
                    s--;
                    d++;
                    while (--off > prev_off) {
                        path.push_back(idx4(qni, tni, off, diag+off));
                        if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);
                    }
                    path.push_back(idx4(qni, tni, off, diag+off+1));
                    if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off+1);
                    diag++;
                    break;
                case PTR_DEL:
                    if (print) printf("deletion\n");
                    prev_off = offs[node_offs.at(prev_node) + d-1];
                    s--;
                    d--;
                    while (--off >= prev_off) {
                        path.push_back(idx4(qni, tni, off, diag+off));
                        if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off, diag+off);
                    }
                    path.push_back(idx4(qni, tni, off+1, diag+off));
                    if (print) printf("path {%d, %d} (%d, %d)\n", qni, tni, off+1, diag+off);
                    diag--;
                    break;
                default:
                    ERROR("Unexpected pointer value at {%d, %d} (%d, %d) "
                            "(prev=(%u, %u), ptr=%u) during WFA backtrack.", 
                            qni, tni, off, diag + off,
                            (ptr & QNODE_MASK) >> (PTR_BITS + NODE_BITS),
                            (ptr & TNODE_MASK) >> PTR_BITS,
                            ptr & PTR_MASK);
                    break;
            }
        }
        off = offs[node_offs.at(idx3(qni, tni, s)) + d];
        ptr = ptrs[node_offs.at(idx3(qni, tni, s)) + d];
    }
    std::reverse(path.begin(), path.end());
    return path;
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

        // init empty pointers/offsets
        std::vector<uint32_t> wfa_ptrs;
        std::vector<int> wfa_offs;
        std::unordered_map<idx3, uint64_t> node_offs;

        int wfa_score = wfa_calc_prec_recall_aln(graph, wfa_ptrs, wfa_offs, node_offs, print);
        bool aligned = wfa_score <= g.max_dist;
        if (aligned) { // alignment succeeded
            std::vector<idx4> path = parse_wfa_path(graph, wfa_score, wfa_ptrs, wfa_offs, node_offs, print);
            calc_prec_recall(graph, path, truth_hi, print);
            done = true;
        }

        // clear pointer and offset matrices as soon as possible
        wfa_ptrs.clear();
        wfa_offs.clear();
        wfa_ptrs.shrink_to_fit();
        wfa_offs.shrink_to_fit();

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
                std::vector< std::vector< std::vector< std::vector<int> > > > retry_offs;
                int dist = gwfa_aln(retry_graph, retry_offs, print);
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
        std::vector<idx4> & path,
        int truth_hap, bool print
        ) {

    idx4 curr = path.back();
    path.pop_back();
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

        idx4 prev = path.back();
        path.pop_back();

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
            && prev.qi+1 == curr.qi;
        bool ref_truth_move = graph->ttypes[curr.tni] == TYPE_REF 
            && prev.tni == curr.tni
            && prev.ti+1 == curr.ti;
        bool sync_point = on_main_diag && (
                (same_submatrix && ref_query_move && ref_truth_move) || diff_submatrix);
        if (sync_point) {
            if (print) printf("potential sync point\n");

            // add sync point
            if (sync_tvars.size() || sync_qvars.size()) {

                // calculate edit distance without variants
                int ref_dist = 0;
                std::vector< std::vector<int> > ed_offs, ed_ptrs;
                int truth_pos = graph->get_truth_pos(curr.tni, curr.ti);
                wf_ed(graph->ref.substr(query_ref_pos, prev_query_ref_pos - query_ref_pos), 
                        graph->truth.substr(truth_pos, prev_truth_pos - truth_pos),
                        ref_dist, ed_offs, ed_ptrs);
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
        std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        std::vector< std::vector< std::vector<int> > > & offs,
        int & s, int x, int o, int e, bool print
        ) {

    // init
    int query_len = query.size();
    int truth_len = truth.size();
    int mat_len = query_len + truth_len - 1;
    bool done = false;
    for (int m = 0; m < MATS; m++) {
        offs[m].push_back(std::vector<int>(mat_len, -2));
        ptrs[m].push_back(std::vector<uint8_t>(mat_len, PTR_NONE));
    }
    s = 0;
    offs[MAT_SUB][s][query_len-1] = -1; // main diag
    ptrs[MAT_SUB][s][query_len-1] = PTR_MAT;

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
                    ptrs[MAT_SUB][s][d] |= (m == MAT_INS) ? PTR_INS : PTR_DEL;
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
            ptrs[m].push_back(std::vector<uint8_t>(mat_len, PTR_NONE));
        }

        for (int d = 0; d < mat_len; d++) {
            int diag = d + 1 - query_len;

            // sub (in SUB)
            if (s-x >= 0 && offs[MAT_SUB][s-x][d] != -2 && 
                            offs[MAT_SUB][s-x][d]+1 < query_len &&
                     diag + offs[MAT_SUB][s-x][d]+1 < truth_len &&
                            offs[MAT_SUB][s-x][d]+1 >= offs[MAT_SUB][s][d]) {
                offs[MAT_SUB][s][d] = offs[MAT_SUB][s-x][d] + 1;
                ptrs[MAT_SUB][s][d] |= PTR_SUB;
                if(print) printf("(S, %d, %d) sub\n", offs[MAT_SUB][s][d], 
                        offs[MAT_SUB][s][d]+diag);
            }

            // open gap (enter DEL)
            if (s-(o+e) >= 0 && d > 0 && 
                           offs[MAT_SUB][s-(o+e)][d-1] != -2 &&
                    diag + offs[MAT_SUB][s-(o+e)][d-1] < truth_len &&
                           offs[MAT_SUB][s-(o+e)][d-1] >= offs[MAT_DEL][s][d]) {
                offs[MAT_DEL][s][d] = offs[MAT_SUB][s-(o+e)][d-1];
                ptrs[MAT_DEL][s][d] |= PTR_SUB;
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
                ptrs[MAT_INS][s][d] |= PTR_SUB;
                if(print) printf("(I, %d, %d) open\n", offs[MAT_INS][s][d], 
                        offs[MAT_INS][s][d]+diag);
            }

            // extend gap (stay DEL)
            if (s-e >= 0 && d > 0 && 
                           offs[MAT_DEL][s-e][d-1] != -2 &&
                    diag + offs[MAT_DEL][s-e][d-1] < truth_len &&
                           offs[MAT_DEL][s-e][d-1] >= offs[MAT_DEL][s][d]) {
                offs[MAT_DEL][s][d] = offs[MAT_DEL][s-e][d-1];
                ptrs[MAT_DEL][s][d] |= PTR_DEL;
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
                ptrs[MAT_INS][s][d] |= PTR_INS;
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


int calc_cig_swg_score(const std::vector<int> & cigar, 
        int sub, int open, int extend) {
    int score = 0;
    int prev = PTR_MAT;
    for (int i = 0; i < int(cigar.size()); i++) {
        switch (cigar[i]) {
            case PTR_MAT:
                i++; // increment by 2
                prev = PTR_MAT;
                break;
            case PTR_SUB:
                score += sub;
                prev = PTR_SUB;
                i++; // increment by 2
                break;
            case PTR_INS:
                score += (prev == PTR_INS ? extend : open + extend);
                prev = PTR_INS;
                break;
            case PTR_DEL:
                score += (prev == PTR_DEL ? extend : open + extend);
                prev = PTR_DEL;
                break;
            default:
                ERROR("Unexpected variant type in calc_cig_swg_score()");
        }
    }
    return score;
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


/* After alignment, backtrack using pointers and save cigar. */
std::vector<int> wf_swg_backtrack(
        const std::string & query,
        const std::string & ref,
        const std::vector< std::vector< std::vector<uint8_t> > > & ptrs, 
        const std::vector< std::vector< std::vector<int> > > & offs, 
        int s, int x, int o, int e,
        bool print) {

    /* if (print) print_wfa_ptrs(query, ref, s, ptrs, offs); */
    /* for (int mi = 0; mi < MATS; mi++) { */
    /*     if(print) printf("\n%s matrix\n", type_strs[mi+1].data()); */
    /*     for (int si = 0; si <= s; si++) { */
    /*         if(print) printf("%d:\n", si); */
    /*         for (int di = 0; di < int(query.size() + ref.size()-1); di++) { */
    /*             if(print) printf("\t%d", offs[mi][si][di]); */
    /*         } */
    /*         if(print) printf("\n"); */
    /*     } */
    /* } */

    // initialize backtrack position and cigar string
    std::vector<int> cigar(query.size() + ref.size(), 0);
    int cigar_ptr = cigar.size()-1;
    idx3 pos(MAT_SUB, query.size()-1, ref.size()-1);

    if (print) printf("\nWF SWG Backtrack:\n");
    while (pos.qi >= 0 || pos.ri >= 0) {

        // debug print
        if (s < 0) ERROR("Negative wf_swg_backtrack() score at (%c, %d, %d)",
                std::string("SID")[pos.mi], pos.qi, pos.ri);

        // init
        int diag = pos.ri - pos.qi;
        int d = query.size()-1 + diag;
        
        if (pos.mi == MAT_SUB) {
            if (ptrs[MAT_SUB][s][d] & PTR_INS) { // end INS freely
                int prev_off = offs[MAT_INS][s][d];
                while (pos.qi > prev_off) { // match
                    if (print) printf("(%c, %d, %d) match ins\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                    if (pos.qi < 0 || pos.ri < 0) ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                }
                pos.mi = MAT_INS;
            } else if (ptrs[MAT_SUB][s][d] & PTR_DEL) { // end DEL freely
                int prev_off = offs[MAT_DEL][s][d];
                while (pos.qi > prev_off) { // match
                    if (print) printf("(%c, %d, %d) match del\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                    if (pos.qi < 0 || pos.ri < 0) ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                }
                pos.mi = MAT_DEL;
            } else if (ptrs[MAT_SUB][s][d] & PTR_SUB) { // sub
                if (s-x < 0) ERROR("Unexpected PTR_SUB");
                int prev_off = offs[MAT_SUB][s-x][d];
                while (pos.qi > prev_off+1) { // match
                    if (print) printf("(%c, %d, %d) match sub\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                    if (pos.qi < 0 || pos.ri < 0) ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                }
                if (print) printf("(%c, %d, %d) sub\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.qi--; cigar[cigar_ptr--] = PTR_SUB;
                pos.ri--; cigar[cigar_ptr--] = PTR_SUB;
                s -= x;
            } else if (ptrs[MAT_SUB][s][d] & PTR_MAT) { // only matches remain
                while (pos.qi >= 0 && pos.ri >= 0) { // match
                    if (print) printf("(%c, %d, %d) match mat\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                }
                if (pos.qi >= 0 || pos.ri >= 0)
                    ERROR("PTR_MAT expected all matches, now at (%c, %d, %d)",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
            } else {
                ERROR("Unexpected pointer '%d' in wf_swg_backtrack() at (%c, %d, %d)",
                    ptrs[MAT_SUB][s][d], std::string("SID")[pos.mi], pos.qi, pos.ri);
            }

        } else if (pos.mi == MAT_INS) {
            if (ptrs[MAT_INS][s][d] & PTR_INS) { // extend INS
                if (print) printf("(%c, %d, %d) ext ins\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.qi--; cigar[cigar_ptr--] = PTR_INS;
                s -= e;
            } else if (ptrs[MAT_INS][s][d] & PTR_SUB) { // start INS
                if (print) printf("(%c, %d, %d) start ins\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.qi--; cigar[cigar_ptr--] = PTR_INS;
                pos.mi = MAT_SUB;
                s -= o + e;
            } else {
                ERROR("Unexpected pointer '%d' in wf_swg_backtrack() at (%c, %d, %d)",
                    ptrs[MAT_INS][s][d], std::string("SID")[pos.mi], pos.qi, pos.ri);
            }

        } else if (pos.mi == MAT_DEL) {
            if (ptrs[MAT_DEL][s][d] & PTR_DEL) { // extend DEL
                if (print) printf("(%c, %d, %d) ext del\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.ri--; cigar[cigar_ptr--] = PTR_DEL;
                s -= e;
            } else if (ptrs[MAT_DEL][s][d] & PTR_SUB) { // start DEL
                if (print) printf("(%c, %d, %d) start del\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.ri--; cigar[cigar_ptr--] = PTR_DEL;
                pos.mi = MAT_SUB;
                s -= o + e;
            } else {
                ERROR("Unexpected pointer '%d' in wf_swg_backtrack() at (%c, %d, %d)",
                    ptrs[MAT_DEL][s][d], std::string("SID")[pos.mi], pos.qi, pos.ri);
            }

        } else {
            ERROR("Unexpected MAT type in wf_swg_backtrack()");
        }
        if (!(pos.qi == -1 && pos.ri == -1) && (pos.qi < 0 || pos.ri < 0)) 
            ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                std::string("SID")[pos.mi], pos.qi, pos.ri);
    }
    return cigar;
}


/******************************************************************************/


void wf_ed(
        const std::string & query, const std::string & truth, int & s,
        std::vector< std::vector<int> > & offs,
        std::vector< std::vector<int> > & ptrs, bool print
        ) {

    // alignment
    int query_len = query.size();
    int truth_len = truth.size();

    // early exit if either string is empty
    if (!query_len) { s = truth_len; return; }
    if (!truth_len) { s = query_len; return; }
    s = 0;

    int mat_len = query_len + truth_len - 1;
    offs.push_back(std::vector<int>(mat_len, -2));
    offs[0][query_len-1] = -1;
    ptrs.push_back(std::vector<int>(mat_len, PTR_NONE));
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
        ptrs.push_back(std::vector<int>(mat_len, PTR_NONE));
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
                ptrs[s][d] |= PTR_SUB;
                if(print) printf("(%d, %d) sub\n", offs[s][d], offs[s][d]+diag);
            }

            // DEL
            if (s-1 >= 0 && d > 0 &&
                           offs[s-1][d-1] != -2 &&
                    diag + offs[s-1][d-1] < truth_len &&
                           offs[s-1][d-1] >= offs[s][d]) {
                offs[s][d] = offs[s-1][d-1];
                ptrs[s][d] |= PTR_DEL;
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
                ptrs[s][d] |= PTR_INS;
                if(print) printf("(%d, %d) ins\n", offs[s][d], offs[s][d]+diag);
            }
        }
    }
}

