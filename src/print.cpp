#include <string>
#include <vector>
#include <cmath>

#include "print.h"
#include "dist.h"
#include "edit.h"

std::string GREEN(int i) { return "\033[32m" + std::to_string(i) + "\033[0m"; }
std::string GREEN(char c) { return "\033[32m" + std::string(1,c) + "\033[0m"; }
std::string GREEN(std::string str) { return "\033[32m" + str + "\033[0m"; }
std::string RED(int i) { return "\033[31m" + std::to_string(i) + "\033[0m"; }
std::string RED(char c) { return "\033[31m" + std::string(1,c) + "\033[0m"; }
std::string RED(std::string str) { return "\033[31m" + str + "\033[0m"; }
std::string BLUE(int i) { return "\033[34m" + std::to_string(i) + "\033[0m"; }
std::string BLUE(char c) { return "\033[34m" + std::string(1,c) + "\033[0m"; }
std::string BLUE(std::string str) { return "\033[34m" + str + "\033[0m"; }
std::string YELLOW(int i) { return "\033[33m" + std::to_string(i) + "\033[0m"; }
std::string YELLOW(char c) { return "\033[33m" + std::string(1,c) + "\033[0m"; }
std::string YELLOW(std::string str) { return "\033[33m" + str + "\033[0m"; }
std::string PURPLE(int i) { return "\033[35m" + std::to_string(i) + "\033[0m"; }
std::string PURPLE(char c) { return "\033[35m" + std::string(1,c) + "\033[0m"; }
std::string PURPLE(std::string str) { return "\033[35m" + str + "\033[0m"; }


/*******************************************************************************/


void print_ref_ptrs(std::vector< std::vector<int> > ptrs) {
    for(int i = 0; i < int(ptrs[0].size()); i++) {
        if (ptrs[FLAGS][i] & PTR_VAR_BEG && ptrs[FLAGS][i] & PTR_VAR_END) {
            printf("%s ", YELLOW(ptrs[PTRS][i]).data());
        } else if (ptrs[FLAGS][i] & PTR_VAR_BEG) {
            printf("%s ", GREEN(ptrs[PTRS][i]).data());
        } else if (ptrs[FLAGS][i] & PTR_VAR_END) {
            printf("%s ", RED(ptrs[PTRS][i]).data());
        } else if (ptrs[FLAGS][i] & PTR_VARIANT) {
            printf("%s ", BLUE(ptrs[PTRS][i]).data());
        } else {
            printf("%d ", ptrs[PTRS][i]);
        }
    }
    printf("\n");
}

/*******************************************************************************/


void print_wfa_ptrs(
        std::vector<std::string> query,
        std::vector<std::string> truth,
        std::vector<int> s,
        std::vector< std::vector< std::vector<int> > > offs, 
        std::vector< std::vector< std::vector<int> > > ptrs) {

    printf("\n=======================================================================\n");
    for(int h = 0; h < 4; h++) { // 4 alignments
        printf("\n%s ALIGNMENT (distance %d)\n", 
                aln_strs[h].data(), s[h]);
        int query_lens = query[h].size();
        int truth_lens = truth[h].size();

        // create array
        std::vector< std::vector<char> > ptr_str;
        for (int i = 0; i < query_lens; i++)
            ptr_str.push_back(std::vector<char>(truth_lens, '.'));

        // modify array with pointers
        int mat_len = query_lens + truth_lens - 1;
        for (int si = 0; si <= s[h]; si++) {
            for(int di = 0; di < mat_len; di++) {
                int diag = di + 1 - query_lens;
                int off = offs[h][si][di];
 
                // check that indices are within bounds
                int query_pos = off;
                int truth_pos = diag + off;
                if (query_pos < 0 || truth_pos < 0) continue;
                if (query_pos > query_lens-1 || 
                        truth_pos > truth_lens-1) continue;
 
                if (si == 0 && diag == 0) { // main diag, no prior edits
                    while (query_pos >= 0)
                        ptr_str[query_pos--][truth_pos--] = '\\';
                } else if (query_pos == 0) { // left edge
                    ptr_str[query_pos][truth_pos] = '-';
                } else if (truth_pos == 0) { // top edge
                    ptr_str[query_pos][truth_pos] = '|';
                } else { // follow diagonal
 
                    int top_off = (di < mat_len-1) ? offs[h][si-1][di+1]+1 : -2;
                    int left_off = (di > 0) ? offs[h][si-1][di-1] : -2;
                    int diag_off = offs[h][si-1][di]+1;
                    while (query_pos > 0 && truth_pos > 0 && 
                            query_pos > top_off && 
                            query_pos > left_off && 
                            query_pos > diag_off) {
                        ptr_str[query_pos--][truth_pos--] = '\\';
                    }
                    // check left/up
                    if (query_pos == diag_off) {
                        ptr_str[query_pos][truth_pos] = 'X';
                    } else if (query_pos == top_off) {
                        ptr_str[query_pos][truth_pos] = '|';
                    } else if (query_pos == left_off) {
                        ptr_str[query_pos][truth_pos] = '-';
                    }
                }
            }
        }
  
        // print array
        for (int i = -1; i < query_lens; i++) {
            for (int j = -1; j < truth_lens; j++) {
                if (i < 0 && j < 0) {
                    printf("  ");
                }
                else if (i < 0) {
                    printf("%c", truth[h][j]);
                } else if (j < 0) {
                    printf("\n%c ", query[h][i]);
                } else {
                    printf("%c", ptr_str[i][j]);
                }
            }
        }
        printf("\n");
    } // 4 alignments
}
           

/*******************************************************************************/


void print_ptrs(std::vector< std::vector<int> > ptrs, 
        std::string alt_str, std::string ref_str) 
{

    // create array
    int alt_len = alt_str.size();
    int ref_len = ref_str.size();
    std::vector< std::vector<char> > ptr_str;
    for (int i = 0; i < alt_len*2; i++)
        ptr_str.push_back(std::vector<char>(ref_len*2, ' '));

    // set arrows
    for (int alt_idx = 0; alt_idx < alt_len; alt_idx++) {
        for(int ref_idx = 0; ref_idx < ref_len; ref_idx++) {
            // states
            if (ptrs[alt_idx][ref_idx] & PTR_SYNC)
                ptr_str[alt_idx*2+1][ref_idx*2+1] = '#';
            else if (ptrs[alt_idx][ref_idx] & PTR_SWP_MAT)
                ptr_str[alt_idx*2+1][ref_idx*2+1] = '+';
            else
                ptr_str[alt_idx*2+1][ref_idx*2+1] = '*';

            // movements
            if (ptrs[alt_idx][ref_idx] & PTR_MAT)
                ptr_str[alt_idx*2][ref_idx*2] = '\\';
            if (ptrs[alt_idx][ref_idx] & PTR_SUB)
                ptr_str[alt_idx*2][ref_idx*2] = 'X';
            if (ptrs[alt_idx][ref_idx] & PTR_DEL)
                ptr_str[alt_idx*2+1][ref_idx*2] = '-';
            if (ptrs[alt_idx][ref_idx] & PTR_INS)
                ptr_str[alt_idx*2][ref_idx*2+1] = '|';
        }
    }

    // print array
    for (int i = -1; i < alt_len*2; i++) {
        for (int j = -1; j < ref_len*2; j++) {
            if (i < 0 && j < 0) {
                printf("\n  ");
            }
            else if (i < 0) {
                if (!(j%2)) printf("%c ", ref_str[j>>1]);
            } else if (j < 0) {
                if (!(i%2)) printf("%c ", alt_str[i>>1]);
                else printf("  ");
            } else {
                switch(ptr_str[i][j]) {
                    case '#':
                    case '+':
                    case '*':
                        if (ptrs[i>>1][j>>1] & PTR_LPATH) {
                            if (ptrs[i>>1][j>>1] & PTR_RPATH) { // both
                                printf("%s", GREEN(ptr_str[i][j]).data());
                            } else { // left
                                printf("%s", YELLOW(ptr_str[i][j]).data());
                            }
                        } else if (ptrs[i>>1][j>>1] & PTR_RPATH) { // right
                            printf("%s", BLUE(ptr_str[i][j]).data());
                        } else { // none
                            printf("%c", ptr_str[i][j]);
                        }
                        break;
                    case 'X':
                    case '\\':
                        if (ptrs[i>>1][j>>1] & PTR_LPATH && 
                                (i>>1) > 0 && (j>>1) > 0 &&
                                ptrs[(i>>1)-1][(j>>1)-1] & PTR_LPATH) {
                            if (ptrs[i>>1][j>>1] & PTR_RPATH &&
                                    ptrs[(i>>1)-1][(j>>1)-1] & PTR_RPATH) { // both
                                printf("%s", GREEN(ptr_str[i][j]).data());
                            } else { // left
                                printf("%s", YELLOW(ptr_str[i][j]).data());
                            }
                        } else if (ptrs[i>>1][j>>1] & PTR_RPATH &&
                                (i>>1) > 0 && (j>>1) > 0 &&
                                ptrs[(i>>1)-1][(j>>1)-1] & PTR_RPATH) { // right
                            printf("%s", BLUE(ptr_str[i][j]).data());
                        } else { // none
                            if (i>>1 == 0 && j>>1 == 0) {
                                    if (ptrs[i>>1][j>>1] & PTR_RPATH &&
                                        ptrs[i>>1][j>>1] & PTR_LPATH)
                                        printf("%s", GREEN(ptr_str[i][j]).data());
                                    else if (ptrs[i>>1][j>>1] & PTR_RPATH)
                                        printf("%s", BLUE(ptr_str[i][j]).data());
                                    else if (ptrs[i>>1][j>>1] & PTR_LPATH)
                                        printf("%s", YELLOW(ptr_str[i][j]).data());
                                    else
                                        printf("%c", ptr_str[i][j]);
                            } else {
                                printf("%c", ptr_str[i][j]);
                            }
                        }
                        break;
                    case '-':
                        if (ptrs[i>>1][j>>1] & PTR_LPATH && 
                                (j>>1) > 0 &&
                                ptrs[i>>1][(j>>1)-1] & PTR_LPATH) {
                            if (ptrs[i>>1][j>>1] & PTR_RPATH &&
                                    ptrs[i>>1][(j>>1)-1] & PTR_RPATH) { // both
                                printf("%s", GREEN('-').data());
                            } else { // left
                                printf("%s", YELLOW('-').data());
                            }
                        } else if (ptrs[i>>1][j>>1] & PTR_RPATH && 
                                (j>>1) > 0 &&
                                ptrs[i>>1][(j>>1)-1] & PTR_RPATH) { // right
                            printf("%s", BLUE('-').data());
                        } else { // none
                            printf("-");
                        }
                        break;
                    case '|':
                        if (ptrs[i>>1][j>>1] & PTR_LPATH && 
                                (i>>1) > 0 &&
                                ptrs[(i>>1)-1][j>>1] & PTR_LPATH) {
                            if (ptrs[i>>1][j>>1] & PTR_RPATH &&
                                    ptrs[(i>>1)-1][j>>1] & PTR_RPATH) { // both
                                printf("%s", GREEN('|').data());
                            } else { // left
                                printf("%s", YELLOW('|').data());
                            }
                        } else if (ptrs[i>>1][j>>1] & PTR_RPATH && 
                                (i>>1) > 0 &&
                                ptrs[(i>>1)-1][j>>1] & PTR_RPATH) { // right
                            printf("%s", BLUE('|').data());
                        } else { // none
                            printf("|");
                        }
                        break;
                    case ' ':
                            printf("%c", ptr_str[i][j]);
                        break;
                }
            }
        }
        printf("\n");
    }
}


/*******************************************************************************/


void write_precision_recall(std::unique_ptr<phaseData> & phasedata_ptr) {

    // init counters; ax0: SUB/INDEL, ax1: TP,FP,FN,PP,PP_FRAC, ax2: QUAL
    int PP_FRAC = 4;
    std::vector< std::vector< std::vector<float> > > query_counts(VARTYPES,
            std::vector< std::vector<float> >(5, 
            std::vector<float>(g.max_qual-g.min_qual+1, 0.0))) ;
    std::vector< std::vector< std::vector<float> > > truth_counts(VARTYPES,
            std::vector< std::vector<float> >(5, 
            std::vector<float>(g.max_qual-g.min_qual+1, 0.0))) ;

    // calculate summary statistics
    for (auto ctg_name : phasedata_ptr->contigs) {

        // add query
        auto & vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->ctg_variants;
        for (int h = 0; h < HAPS; h++) {
            for (int i = 0; i < vars[QUERY][h]->n; i++) {
                float q = vars[QUERY][h]->callq[i];
                int t = 0;
                if (vars[QUERY][h]->types[i] == TYPE_SUB) {
                    t = VARTYPE_SNP;
                } else if (vars[QUERY][h]->types[i] == TYPE_INS ||
                        vars[QUERY][h]->types[i] == TYPE_DEL) {
                    t = VARTYPE_INDEL;
                } else {
                    ERROR("Unexpected variant type (%d) in write_precision_recall()", 
                            vars[QUERY][h]->types[i]);
                }
                if (vars[QUERY][h]->errtypes[i] == ERRTYPE_UN) {
                    WARN("Unknown error type at hap %d variant %d", h+1, i);
                    continue;
                }
                for (int qual = g.min_qual; qual <= q; qual++) {
                    query_counts[t][ vars[QUERY][h]->errtypes[i] ][qual-g.min_qual]++;
                    if (vars[QUERY][h]->errtypes[i] == ERRTYPE_PP) {
                        query_counts[t][PP_FRAC][qual-g.min_qual] += vars[QUERY][h]->credit[i];
                    }
                }
            }
        }

        // add truth
        for (int h = 0; h < HAPS; h++) {
            for (int i = 0; i < vars[TRUTH][h]->n; i++) {
                float q = vars[TRUTH][h]->callq[i];
                int t = 0;
                if (vars[TRUTH][h]->types[i] == TYPE_SUB) {
                    t = VARTYPE_SNP;
                } else if (vars[TRUTH][h]->types[i] == TYPE_INS ||
                        vars[TRUTH][h]->types[i] == TYPE_DEL) {
                    t = VARTYPE_INDEL;
                } else {
                    ERROR("Unexpected variant type (%d) in write_precision_recall()", 
                            vars[TRUTH][h]->types[i]);
                }
                for (int qual = g.min_qual; qual <= q; qual++) {
                    truth_counts[t][ vars[TRUTH][h]->errtypes[i] ][qual-g.min_qual]++;
                    if (vars[TRUTH][h]->errtypes[i] == ERRTYPE_PP) {
                        truth_counts[t][PP_FRAC][qual-g.min_qual] += vars[TRUTH][h]->credit[i];
                    }
                }
            }
        }
    }

    // write results
    std::string out_pr_fn = g.out_prefix + "precision-recall.tsv";
    INFO("  Printing precision-recall results to '%s'", out_pr_fn.data());
    FILE* out_pr = fopen(out_pr_fn.data(), "w");
    fprintf(out_pr, "VAR_TYPE\tMIN_QUAL\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE\t"
            "TRUTH_TOTAL\tTRUTH_TP\tTRUTH_PP\tTRUTH_FN\tQUERY_TOTAL\tQUERY_TP\tQUERY_PP\tQUERY_FP\n");
    std::vector<float> max_f1_score = {0, 0};
    std::vector<int> max_f1_qual = {0, 0};
    for (int type = 0; type < VARTYPES; type++) {

        // only sweeping query qualities; always consider all truth variants
        for (int qual = g.min_qual; qual <= g.max_qual; qual++) {
            int qidx = qual - g.min_qual;

            int query_tot = query_counts[type][ERRTYPE_TP][qidx] +
                        query_counts[type][ERRTYPE_FP][qidx] +
                        query_counts[type][ERRTYPE_PP][qidx];
            float query_tp_f = query_counts[type][ERRTYPE_TP][qidx] +
                        query_counts[type][PP_FRAC][qidx];
            float query_fp_f = query_tot - query_tp_f;

            int truth_tot = truth_counts[type][ERRTYPE_TP][0] +
                        truth_counts[type][ERRTYPE_PP][0] +
                        truth_counts[type][ERRTYPE_FN][0];
            float truth_tp_f = truth_counts[type][ERRTYPE_TP][qidx] +
                         truth_counts[type][PP_FRAC][qidx];
            if (truth_tot == 0) break;
            if (truth_tp_f + query_fp_f == 0) break;

            // ignore PP, this is only for summary output, not calculations
            int truth_fn = truth_counts[type][ERRTYPE_FN][0] +
                         truth_counts[type][ERRTYPE_TP][0] -
                         truth_counts[type][ERRTYPE_TP][qidx];

            float precision = truth_tp_f / (truth_tp_f + query_fp_f);
            float recall = truth_tp_f / truth_tot;
            float f1_score = 2*precision*recall / (precision + recall);
            if (f1_score > max_f1_score[type]) {
                max_f1_score[type] = f1_score;
                max_f1_qual[type] = qual;
            }

            fprintf(out_pr, 
                    "%s\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    vartype_strs[type].data(),
                    qual,
                    precision,
                    recall, 
                    f1_score,
                    qscore(1-f1_score),
                    truth_tot,
                    int(truth_counts[type][ERRTYPE_TP][qidx]),
                    int(truth_counts[type][ERRTYPE_PP][qidx]),
                    truth_fn,
                    query_tot, 
                    int(query_counts[type][ERRTYPE_TP][qidx]),
                    int(query_counts[type][ERRTYPE_PP][qidx]),
                    int(query_counts[type][ERRTYPE_FP][qidx])
           );
        }
    }
    fclose(out_pr);

    // print summary output
    std::string out_pr_summ_fn = g.out_prefix + "precision-recall-summary.tsv";
    INFO("  Printing precision-recall summary to '%s'", out_pr_summ_fn.data());
    FILE* out_pr_summ = fopen(out_pr_summ_fn.data(), "w");
    fprintf(out_pr_summ, "VAR_TYPE\tMIN_QUAL\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\t\tRECALL\t\tF1_SCORE\tF1_QSCORE\n");
    for (int type = 0; type < VARTYPES; type++) {
        std::vector<int> quals = {g.min_qual, max_f1_qual[type]};
        INFO(" ");
        INFO("TYPE\tMIN_QUAL\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\t\tRECALL\t\tF1_SCORE\tF1_QSCORE");

        for (int qual : quals) {
            // redo calculations for these two
            int qidx = qual - g.min_qual;

            int query_tot = query_counts[type][ERRTYPE_TP][qidx] +
                        query_counts[type][ERRTYPE_FP][qidx] +
                        query_counts[type][ERRTYPE_PP][qidx];
            float query_tp_f = query_counts[type][ERRTYPE_TP][qidx] +
                         query_counts[type][PP_FRAC][qidx];
            float query_fp_f = query_tot - query_tp_f;
            if (query_tot == 0) WARN("No QUERY variant calls.");

            int truth_tot = truth_counts[type][ERRTYPE_TP][0] +
                        truth_counts[type][ERRTYPE_PP][0] +
                        truth_counts[type][ERRTYPE_FN][0];
            int truth_fn = truth_counts[type][ERRTYPE_FN][0] +
                             truth_counts[type][ERRTYPE_TP][0] -
                             truth_counts[type][ERRTYPE_TP][qidx];
            float truth_tp_f = truth_counts[type][ERRTYPE_TP][qidx] +
                             truth_counts[type][PP_FRAC][qidx];
            if (truth_tot == 0) WARN("No TRUTH variant calls.");

            float precision = truth_tp_f + query_fp_f == 0 ? 
                1.0f : truth_tp_f / (truth_tp_f + query_fp_f);
            float recall = truth_tot == 0 ? 1.0f : truth_tp_f / truth_tot;
            float f1_score = 2*precision*recall / (precision + recall);

            // print summary
            INFO("%s\tQ >= %d\t\t%-16d%-16d%-16d%-16d%f\t%f\t%f\t%f",
                vartype_strs[type].data(),
                qual,
                int(truth_counts[type][ERRTYPE_TP][qidx]),
                int(query_counts[type][ERRTYPE_TP][qidx]),
                truth_fn,
                int(query_counts[type][ERRTYPE_FP][qidx]),
                precision,
                recall, 
                f1_score,
                qscore(1-f1_score)
            );
            fprintf(out_pr_summ,
               "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",
                vartype_strs[type].data(),
                qual,
                int(truth_counts[type][ERRTYPE_TP][qidx]),
                int(query_counts[type][ERRTYPE_TP][qidx]),
                truth_fn,
                int(query_counts[type][ERRTYPE_FP][qidx]),
                precision,
                recall, 
                f1_score,
                qscore(1-f1_score)
            );
        }
    }
    INFO(" ");
    fclose(out_pr_summ);
}


/*******************************************************************************/


void write_distance(const editData & edits) {

    // log all distance results
    std::string dist_fn = g.out_prefix + "distance.tsv";
    FILE* out_dists = fopen(dist_fn.data(), "w");
    fprintf(out_dists, "MIN_QUAL\tSUB_DE\tINS_DE\tDEL_DE\tSUB_ED\tINS_ED\tDEL_ED\t"
            "DISTINCT_EDITS\tEDIT_DIST\tALN_SCORE\tALN_QSCORE\n");
    INFO("  Printing distance results to '%s'", dist_fn.data());

    // get original scores / distance (above g.max_qual, no vars applied)
    std::vector<int> orig_edit_dists(TYPES, 0);
    std::vector<int> orig_distinct_edits(TYPES, 0);
    int orig_score = edits.get_score(g.max_qual+1);
    for (int type = 0; type < TYPES; type++) {
        orig_edit_dists[type] = edits.get_ed(g.max_qual+1, type);
        orig_distinct_edits[type] = edits.get_de(g.max_qual+1, type);
    }

    std::vector<double> best_score(TYPES, std::numeric_limits<double>::max());
    std::vector<int> best_qual(TYPES, 0);
    for (int q = g.min_qual; q <= g.max_qual+1; q++) { // all qualities

        // get ED/DE for each Q threshold, for each type
        std::vector<int> edit_dists(TYPES, 0);
        std::vector<int> distinct_edits(TYPES, 0);
        for (int type = 0; type < TYPES; type++) {
            edit_dists[type] = edits.get_ed(q, type);
            distinct_edits[type] = edits.get_de(q, type);

            // save best Q threshold so far
            double score = double(edit_dists[type]) * distinct_edits[type];
            if (score < best_score[type]) {
                best_score[type] = score;
                best_qual[type] = q;
            }
        }

        fprintf(out_dists, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", 
                q, distinct_edits[TYPE_SUB],
                distinct_edits[TYPE_INS], distinct_edits[TYPE_DEL],
                edit_dists[TYPE_SUB], edit_dists[TYPE_INS], edit_dists[TYPE_DEL],
                distinct_edits[TYPE_ALL], edit_dists[TYPE_ALL],
                edits.get_score(q), qscore(double(edits.get_score(q))/orig_score));
    }
    fclose(out_dists);

    // summarize distance results
    std::string dist_summ_fn = g.out_prefix + "distance-summary.tsv";
    FILE* dists_summ = fopen(dist_summ_fn.data(), "w");
    INFO("  Printing distance summary to '%s'", dist_summ_fn.data());
    fprintf(dists_summ, "VAR_TYPE\tMIN_QUAL\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE\tALN_QSCORE\n");
    for (int type = 0; type < TYPES; type++) {

        // skip INS/DEL individually unless higher print verbosity
        if ((type == TYPE_INS || type == TYPE_DEL) && g.verbosity == 0)
            continue;

        INFO(" ");
        if (type == TYPE_ALL) {
            INFO("TYPE\tMIN_QUAL\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE\tALN_QSCORE");
        } else {
            INFO("TYPE\tMIN_QUAL\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE");
        }
        std::vector<int> quals = {g.min_qual, best_qual[type], g.max_qual+1};
        for (auto q : quals) {

            // fill out ED/DE for selected quals
            std::vector<int> edit_dists(TYPES, 0);
            std::vector<int> distinct_edits(TYPES, 0);
            for (int type = 0; type < TYPES; type++) {
                edit_dists[type] = edits.get_ed(q, type);
                distinct_edits[type] = edits.get_de(q, type);
            }

            float ed_qscore = qscore(double(edit_dists[type]) / orig_edit_dists[type]);
            float de_qscore = qscore(double(distinct_edits[type]) / orig_distinct_edits[type]);
            float all_qscore = type == TYPE_ALL ? qscore(double(edits.get_score(q)) / orig_score) : 0;

            // print summary
            fprintf(dists_summ, "%s\t%d\t%d\t%d\t%f\t%f\t%f\n", type_strs2[type].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, all_qscore);
            if (type == TYPE_ALL) {
                INFO("%s\tQ >= %d\t\t%-16d%-16d%f\t%f\t%f", type_strs2[type].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, all_qscore);
            } else {
                INFO("%s\tQ >= %d\t\t%-16d%-16d%f\t%f", type_strs2[type].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore);
            }
        }
    }
    INFO(" ");
    fclose(dists_summ);
}


/*******************************************************************************/


void write_results(
        std::unique_ptr<phaseData> & phasedata_ptr, 
        const editData & edits) {
    INFO(" ");
    INFO("Writing results");

    // print summary (precision/recall) information
    write_precision_recall(phasedata_ptr);

    // print distance information
    write_distance(edits);

    // print edit information
    std::string edit_fn = g.out_prefix + "edits.tsv";
    FILE* out_edits = fopen(edit_fn.data(), "w");
    fprintf(out_edits, "CONTIG\tSTART\tHAP\tTYPE\tSIZE\tSUPERCLUSTER\tQUAL\n");
    INFO("  Printing edit results to '%s'", edit_fn.data());
    for (int i = 0; i < edits.n; i++) {
        fprintf(out_edits, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n", 
                edits.ctgs[i].data(), edits.poss[i], edits.haps[i],
                type_strs[edits.types[i]].data(), edits.lens[i],
                edits.superclusters[i], edits.quals[i]);
    }
    fclose(out_edits);

    // print phasing information
    std::string out_phasings_fn = g.out_prefix + "phase-blocks.tsv";
    FILE* out_phasings = fopen(out_phasings_fn.data(), "w");
    INFO("  Printing phasing results to '%s'", out_phasings_fn.data());
    fprintf(out_phasings, "CONTIG\tSTART\tSTOP\tSIZE\tSUPERCLUSTERS\n");
    for (auto ctg_name : phasedata_ptr->contigs) {
        ctgPhasings ctg_phasings = phasedata_ptr->ctg_phasings[ctg_name];
        std::shared_ptr<ctgSuperclusters> ctg_superclusters = ctg_phasings.ctg_superclusters;
        for (int i = 0; i <= ctg_phasings.nswitches && ctg_phasings.nswitches > 0; i++) {
            int beg_idx = ctg_phasings.phase_blocks[i];
            int end_idx = ctg_phasings.phase_blocks[i+1]-1;
            int beg = ctg_superclusters->begs[beg_idx];
            int end = ctg_superclusters->ends[end_idx];
            fprintf(out_phasings, "%s\t%d\t%d\t%d\t%d\n", 
                    ctg_name.data(), beg, end, end-beg, end_idx-beg_idx+1);
        }
    }
    fclose(out_phasings);

    // print clustering information
    std::string out_clusterings_fn = g.out_prefix + "superclusters.tsv";
    FILE* out_clusterings = fopen(out_clusterings_fn.data(), "w");
    INFO("  Printing superclustering results to '%s'", out_clusterings_fn.data());
    fprintf(out_clusterings, "CONTIG\tSTART\tSTOP\tSIZE\tQUERY1_VARS\tQUERY2_VARS"
            "\tTRUTH1_VARS\tTRUTH2_VARS\tORIG_ED\tSWAP_ED\tPHASE\tPHASE_BLOCK\n");
    for (auto ctg_name : phasedata_ptr->contigs) {
        auto ctg_phasings = phasedata_ptr->ctg_phasings[ctg_name];
        std::shared_ptr<ctgSuperclusters> ctg_supclusts = ctg_phasings.ctg_superclusters;
        int phase_block_idx = 0;
        for (int i = 0; i < ctg_supclusts->n; i++) {

            // we've entered the next phase block
            if (i >= ctg_phasings.phase_blocks[phase_block_idx+1]) {
                phase_block_idx++;
            }

            // count query vars, allowing empty haps
            int query1_vars = ctg_supclusts->ctg_variants[QUERY][HAP1]->clusters.size() ?
                ctg_supclusts->ctg_variants[QUERY][HAP1]->clusters[
                        ctg_supclusts->superclusters[QUERY][HAP1][i+1]] -
                    ctg_supclusts->ctg_variants[QUERY][HAP1]->clusters[
                        ctg_supclusts->superclusters[QUERY][HAP1][i]] : 0;
            int query2_vars = ctg_supclusts->ctg_variants[QUERY][HAP2]->clusters.size() ?
                ctg_supclusts->ctg_variants[QUERY][HAP2]->clusters[
                        ctg_supclusts->superclusters[QUERY][HAP2][i+1]] -
                    ctg_supclusts->ctg_variants[QUERY][HAP2]->clusters[
                        ctg_supclusts->superclusters[QUERY][HAP2][i]] : 0;
            int truth1_vars = ctg_supclusts->ctg_variants[TRUTH][HAP1]->clusters.size() ?
                ctg_supclusts->ctg_variants[TRUTH][HAP1]->clusters[
                        ctg_supclusts->superclusters[TRUTH][HAP1][i+1]] -
                    ctg_supclusts->ctg_variants[TRUTH][HAP1]->clusters[
                        ctg_supclusts->superclusters[TRUTH][HAP1][i]] : 0;
            int truth2_vars = ctg_supclusts->ctg_variants[TRUTH][HAP2]->clusters.size() ?
                ctg_supclusts->ctg_variants[TRUTH][HAP2]->clusters[
                        ctg_supclusts->superclusters[TRUTH][HAP2][i+1]] -
                    ctg_supclusts->ctg_variants[TRUTH][HAP2]->clusters[
                        ctg_supclusts->superclusters[TRUTH][HAP2][i]] : 0;

            // print data
            fprintf(out_clusterings, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\n", 
                ctg_name.data(), 
                ctg_supclusts->begs[i],
                ctg_supclusts->ends[i],
                ctg_supclusts->ends[i] - ctg_supclusts->begs[i],
                query1_vars, query2_vars, truth1_vars, truth2_vars,
                ctg_supclusts->orig_phase_dist[i],
                ctg_supclusts->swap_phase_dist[i],
                phase_strs[ctg_supclusts->phase[i]].data(),
                phase_block_idx
           );
        }
    }
    fclose(out_clusterings);

    // print query variant information
    std::string out_query_fn = g.out_prefix + "query.tsv";
    INFO("  Printing call variant results to '%s'", out_query_fn.data());
    FILE* out_query = fopen(out_query_fn.data(), "w");
    fprintf(out_query, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERR_TYPE"
            "\tCREDIT\tCLUSTER\tSUPERCLUSTER\tLOCATION\n");
    for (auto ctg_name : phasedata_ptr->contigs) {

        // set pointers to variants and superclusters
        std::shared_ptr<ctgVariants> query1_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->ctg_variants[QUERY][HAP1];
        std::shared_ptr<ctgVariants> query2_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->ctg_variants[QUERY][HAP2];
        std::shared_ptr<ctgSuperclusters> ctg_supclusts = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters;

        int supercluster_idx = 0;
        int var1_idx = 0;
        int var2_idx = 0;
        int cluster1_idx = 0;
        int cluster2_idx = 0;
        while (var1_idx < query1_vars->n || var2_idx < query2_vars->n) {
            if (var2_idx >= query2_vars->n || // only query1 has remaining vars
                    (var1_idx < query1_vars->n && query1_vars->poss[var1_idx] < query2_vars->poss[var2_idx])) { // query1 var next

                // we've entered the next supercluster
                if (cluster1_idx+1 >= int(query1_vars->clusters.size())) 
                    ERROR("Out of bounds cluster during write_results(): query1")
                if (query1_vars->clusters[cluster1_idx+1] <= var1_idx) cluster1_idx++;
                if (supercluster_idx >= int(ctg_supclusts->begs.size())) 
                    ERROR("Out of bounds supercluster during write_results(): query1")
                if (query1_vars->poss[var1_idx] >= ctg_supclusts->ends[supercluster_idx])
                    supercluster_idx++;
                fprintf(out_query, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        query1_vars->poss[var1_idx],
                        query1_vars->haps[var1_idx],
                        query1_vars->refs[var1_idx].data(),
                        query1_vars->alts[var1_idx].data(),
                        query1_vars->var_quals[var1_idx],
                        type_strs[query1_vars->types[var1_idx]].data(),
                        error_strs[query1_vars->errtypes[var1_idx]].data(),
                        query1_vars->credit[var1_idx],
                        cluster1_idx,
                        supercluster_idx,
                        region_strs[query1_vars->locs[var1_idx]].data()
                       );
                var1_idx++;
            } else { // process query2 var
                if (cluster2_idx+1 >= int(query2_vars->clusters.size()))
                    ERROR("Out of bounds cluster during write_results(): query2")
                if (query2_vars->clusters[cluster2_idx+1] <= var2_idx) cluster2_idx++;
                if (supercluster_idx >= int(ctg_supclusts->begs.size())) 
                    ERROR("Out of bounds supercluster during write_results(): query2")
                if (query2_vars->poss[var2_idx] >= ctg_supclusts->ends[supercluster_idx])
                    supercluster_idx++;
                fprintf(out_query, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        query2_vars->poss[var2_idx],
                        query2_vars->haps[var2_idx],
                        query2_vars->refs[var2_idx].data(),
                        query2_vars->alts[var2_idx].data(),
                        query2_vars->var_quals[var2_idx],
                        type_strs[query2_vars->types[var2_idx]].data(),
                        error_strs[query2_vars->errtypes[var2_idx]].data(),
                        query2_vars->credit[var2_idx],
                        cluster2_idx,
                        supercluster_idx,
                        region_strs[query2_vars->locs[var2_idx]].data()
                       );
                var2_idx++;
            }
        }
    }
    fclose(out_query);
    
    //print truth variant information
    std::string out_truth_fn = g.out_prefix + "truth.tsv";
    INFO("  Printing truth variant results to '%s'", out_truth_fn.data());
    FILE* out_truth = fopen(out_truth_fn.data(), "w");
    fprintf(out_truth, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERRTYPE\tCREDIT\tCLUSTER\tSUPERCLUSTER\tLOCATION\n");
    for (auto ctg_name : phasedata_ptr->contigs) {

        // set pointers to variants and superclusters
        std::shared_ptr<ctgVariants> truth1_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->ctg_variants[TRUTH][HAP1];
        std::shared_ptr<ctgVariants> truth2_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->ctg_variants[TRUTH][HAP2];
        std::shared_ptr<ctgSuperclusters> ctg_supclusts = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters;

        int var1_idx = 0;
        int var2_idx = 0;
        int cluster1_idx = 0;
        int cluster2_idx = 0;
        int supercluster_idx = 0;
        while (var1_idx < truth1_vars->n || var2_idx < truth2_vars->n) {
            if (var2_idx >= truth2_vars->n || // only truth1 has remaining vars
                    (var1_idx < truth1_vars->n && truth1_vars->poss[var1_idx] < truth2_vars->poss[var2_idx])) { // truth1 var next
                if (cluster1_idx+1 >= int(truth1_vars->clusters.size()))
                    ERROR("Out of bounds cluster during write_results(): truth1")
                if (truth1_vars->clusters[cluster1_idx+1] <= var1_idx) cluster1_idx++;
                if (supercluster_idx >= int(ctg_supclusts->begs.size())) 
                    ERROR("Out of bounds supercluster during write_results(): truth1")
                if (truth1_vars->poss[var1_idx] >= ctg_supclusts->ends[supercluster_idx])
                    supercluster_idx++;
                fprintf(out_truth, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        truth1_vars->poss[var1_idx],
                        truth1_vars->haps[var1_idx],
                        truth1_vars->refs[var1_idx].data(),
                        truth1_vars->alts[var1_idx].data(),
                        truth1_vars->var_quals[var1_idx],
                        type_strs[truth1_vars->types[var1_idx]].data(),
                        error_strs[truth1_vars->errtypes[var1_idx]].data(),
                        truth1_vars->credit[var1_idx],
                        cluster1_idx,
                        supercluster_idx,
                        region_strs[truth1_vars->locs[var1_idx]].data()
                       );
                var1_idx++;
            } else { // process truth2 var
                if (cluster2_idx+1 >= int(truth2_vars->clusters.size())) {
                    ERROR("Out of bounds cluster during write_results(): truth2")
                }
                if (truth2_vars->clusters[cluster2_idx+1] <= var2_idx) cluster2_idx++;
                if (supercluster_idx >= int(ctg_supclusts->begs.size())) 
                    ERROR("Out of bounds supercluster during write_results(): truth2")
                if (truth2_vars->poss[var2_idx] >= ctg_supclusts->ends[supercluster_idx])
                    supercluster_idx++;
                fprintf(out_truth, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        truth2_vars->poss[var2_idx],
                        truth2_vars->haps[var2_idx],
                        truth2_vars->refs[var2_idx].data(),
                        truth2_vars->alts[var2_idx].data(),
                        truth2_vars->var_quals[var2_idx],
                        type_strs[truth2_vars->types[var2_idx]].data(),
                        error_strs[truth2_vars->errtypes[var2_idx]].data(),
                        truth2_vars->credit[var2_idx],
                        cluster2_idx,
                        supercluster_idx,
                        region_strs[truth2_vars->locs[var2_idx]].data()
                       );
                var2_idx++;
            }
        }
    }
    fclose(out_truth);

    return;
}


/*******************************************************************************/


void print_cigar(std::vector<int> cigar) {
    for (int i = 0; i < int(cigar.size()); i++) {
        switch (cigar[i]) {
        case PTR_MAT:
            i++;
            printf("=");
            break;
        case PTR_SUB:
            i++;
            printf("X");
            break;
        case PTR_INS:
            printf("I");
            break;
        case PTR_DEL:
            printf("D");
            break;
        default:
            printf("?");
            break;
            /* ERROR("Unexpected CIGAR type in print_cigar(): %d", cigar[i]); */
        }
    }
    printf("\n");
}
