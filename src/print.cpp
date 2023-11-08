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

inline std::string b2s(bool b) { return b ? "true" : "false"; }

void write_params() {

    // concatenate filters into string
    std::string filters_str = g.filters.size() ? g.filters[0] : "";
    for (int i = 1; i < int(g.filters.size()); i++)
        filters_str += "," + g.filters[i];

    // write all params to file
    std::string out_params_fn = g.out_prefix + "parameters.txt";
    FILE* out_params = fopen(out_params_fn.data(), "w");
    fprintf(out_params, 
        "program = '%s'\nversion = '%s'\nout_prefix = '%s'\ncommand = '%s'\nreference_fasta = '%s'\n"
        "query_vcf = '%s'\ntruth_vcf = '%s'\nbed_file = '%s'\nwrite_outputs = %s\nfilters = '%s'\n"
        "min_var_qual = %d\nmax_var_qual = %d\nmin_var_size = %d\nmax_var_size = %d\n"
        "phase_threshold = %f\ncredit_threshold = %f\nrealign_truth = %s\nrealign_query = %s\n"
        "realign_only = %s\nsimple_cluster = %s\ncluster_min_gap = %d\n"
        "reach_min_gap = %d\nmax_cluster_itrs = %d\nmax_threads = %d\nmax_ram = %f\n"
        "sub = %d\nopen = %d\nextend = %d\neval_sub = %d\neval_open = %d\neval_extend = %d\n",
        g.PROGRAM.data(), g.VERSION.data(), g.out_prefix.data(), g.cmd.data(), g.ref_fasta_fn.data(), 
        g.query_vcf_fn.data(), g.truth_vcf_fn.data(), g.bed_fn.data(), b2s(g.write).data(), 
        filters_str.data(), g.min_qual, g.max_qual, g.min_size, g.max_size,
        g.phase_threshold, g.credit_threshold, b2s(g.realign_truth).data(), b2s(g.realign_query).data(),
        b2s(g.realign_only).data(), b2s(g.simple_cluster).data(), g.cluster_min_gap,
        g.reach_min_gap, g.max_cluster_itrs, g.max_threads, g.max_ram,
        g.sub, g.open, g.extend, g.eval_sub, g.eval_open, g.eval_extend);
    fclose(out_params);
}


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
        const std::string & query,
        const std::string & truth,
        int s,
        const std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        const std::vector< std::vector< std::vector<int> > > & offs) {

    for (int m = 0; m < MATS; m++) {
        int query_len = query.size();
        int truth_len = truth.size();

        // create array
        std::vector< std::vector<char> > ptr_str;
        for (int i = 0; i < query_len; i++)
            ptr_str.push_back(std::vector<char>(truth_len, '.'));

        // modify array with pointers
        int mat_len = query_len + truth_len - 1;
        for (int si = s; si >= 0; si--) { // overwrite with lower scored ptrs
            for(int di = 0; di < mat_len; di++) {
                int diag = di + 1 - query_len;
                int off = offs[m][si][di];
 
                // check that indices are within bounds
                int query_pos = off;
                int truth_pos = diag + off;
                if (query_pos < 0 || truth_pos < 0) continue;
                if (query_pos > query_len-1 || 
                        truth_pos > truth_len-1) continue;

                if (m == MAT_SUB) {
                    switch (ptrs[m][si][di]) {
                        case PTR_SUB:
                            ptr_str[query_pos][truth_pos] = 'X';
                            break;
                        case PTR_INS:
                            ptr_str[query_pos][truth_pos] = 'i';
                            break;
                        case PTR_DEL:
                            ptr_str[query_pos][truth_pos] = 'd';
                            break;
                        case PTR_MAT:
                            ptr_str[query_pos][truth_pos] = '\\';
                            break;
                    }

                } else if (m == MAT_INS) {
                    switch (ptrs[m][si][di]) {
                        case PTR_SUB:
                            ptr_str[query_pos][truth_pos] = 's';
                            break;
                        case PTR_INS:
                            ptr_str[query_pos][truth_pos] = '|';
                            break;
                        case PTR_DEL:
                            ptr_str[query_pos][truth_pos] = '?';
                            break;
                    }

                } else if (m == MAT_DEL) {
                    switch (ptrs[m][si][di]) {
                        case PTR_SUB:
                            ptr_str[query_pos][truth_pos] = 's';
                            break;
                        case PTR_INS:
                            ptr_str[query_pos][truth_pos] = '?';
                            break;
                        case PTR_DEL:
                            ptr_str[query_pos][truth_pos] = '-';
                            break;
                    }
                }
            }
        }
  
        // print array
        printf("\n%s matrix:\n", type_strs[m+1].data());
        for (int i = -1; i < query_len; i++) {
            for (int j = -1; j < truth_len; j++) {
                if (i < 0 && j < 0) {
                    printf("  ");
                }
                else if (i < 0) {
                    printf("%c", truth[j]);
                } else if (j < 0) {
                    printf("\n%c ", query[i]);
                } else {
                    printf("%c", ptr_str[i][j]);
                }
            }
        }
        printf("\n");
    }
}
           
/*******************************************************************************/

void print_ptrs(const std::vector< std::vector<uint8_t> > & ptrs, 
        const std::string & alt_str, const std::string & ref_str) 
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


void write_precision_recall(std::unique_ptr<phaseblockData> & phasedata_ptr) {

    // for each class, store variant counts above each quality threshold
    // init counters; ax0: SUB/INDEL, ax1: TP,FP,FN ax2: QUAL
    std::vector< std::vector< std::vector<float> > > query_counts(VARTYPES,
            std::vector< std::vector<float> >(3, 
            std::vector<float>(g.max_qual-g.min_qual+1, 0.0))) ;
    std::vector< std::vector< std::vector<float> > > truth_counts(VARTYPES,
            std::vector< std::vector<float> >(3, 
            std::vector<float>(g.max_qual-g.min_qual+1, 0.0))) ;

    // calculate summary statistics
    for (std::string ctg : phasedata_ptr->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
        std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
        auto & ctg_vars = ctg_scs->ctg_variants;

        // add query
        if(ctg_scs->n)
        for (int h = 0; h < HAPS; h++) {
            if (ctg_vars[QUERY][h]->n)
            for (int pbi = 0; pbi < ctg_pbs->n; pbi++) {
                for (int sci = ctg_pbs->phase_blocks[pbi]; 
                        sci < ctg_pbs->phase_blocks[pbi+1]; sci++) {

                    // get swap info
                    bool swap;
                    switch (ctg_scs->phase[sci]) {
                        case PHASE_ORIG: swap = false; break;
                        case PHASE_SWAP: swap = true; break;
                        default: swap = ctg_pbs->block_states[pbi]; break;
                    }

                    for (int i = ctg_vars[QUERY][h]->clusters[
                            ctg_scs->superclusters[QUERY][h][sci] ]; 
                            i < ctg_vars[QUERY][h]->clusters[
                            ctg_scs->superclusters[QUERY][h][sci+1] ]; i++) {

                        float q = ctg_vars[QUERY][h]->callq[swap][i];
                        int t = 0;
                        if (ctg_vars[QUERY][h]->types[i] == TYPE_SUB) {
                            t = VARTYPE_SNP;
                        } else if (ctg_vars[QUERY][h]->types[i] == TYPE_INS ||
                                ctg_vars[QUERY][h]->types[i] == TYPE_DEL) {
                            t = VARTYPE_INDEL;
                        } else {
                            ERROR("Unexpected variant type (%d) in write_precision_recall()", 
                                    ctg_vars[QUERY][h]->types[i]);
                        }
                        if (ctg_vars[QUERY][h]->errtypes[swap][i] == ERRTYPE_UN) {
                            WARN("Unknown error type at QUERY %s:%d", ctg.data(), ctg_vars[QUERY][h]->poss[i]);
                            continue;
                        }
                        for (int qual = g.min_qual; qual <= q; qual++) {
                            query_counts[t][ ctg_vars[QUERY][h]->errtypes[swap][i] ][qual-g.min_qual]++;
                        }
                    }
                }
            }
        }

        // add truth
        if (ctg_scs->n)
        for (int h = 0; h < HAPS; h++) {
            if (ctg_vars[TRUTH][h]->n)
            for (int pbi = 0; pbi < ctg_pbs->n; pbi++) {
                for (int sci = ctg_pbs->phase_blocks[pbi]; 
                        sci < ctg_pbs->phase_blocks[pbi+1]; sci++) {

                    // get swap info
                    bool swap;
                    switch (ctg_scs->phase[sci]) {
                        case PHASE_ORIG: swap = false; break;
                        case PHASE_SWAP: swap = true; break;
                        default: swap = ctg_pbs->block_states[pbi]; break;
                    }

                    for (int i = ctg_vars[TRUTH][h]->clusters[
                            ctg_scs->superclusters[TRUTH][h][sci] ]; 
                            i < ctg_vars[TRUTH][h]->clusters[
                            ctg_scs->superclusters[TRUTH][h][sci+1] ]; i++) {

                        float q = ctg_vars[TRUTH][h]->callq[swap][i];
                        int t = 0;
                        if (ctg_vars[TRUTH][h]->types[i] == TYPE_SUB) {
                            t = VARTYPE_SNP;
                        } else if (ctg_vars[TRUTH][h]->types[i] == TYPE_INS ||
                                ctg_vars[TRUTH][h]->types[i] == TYPE_DEL) {
                            t = VARTYPE_INDEL;
                        } else {
                            ERROR("Unexpected variant type (%d) in write_precision_recall()", 
                                    ctg_vars[TRUTH][h]->types[i]);
                        }
                        if (ctg_vars[TRUTH][h]->errtypes[swap][i] == ERRTYPE_UN) {
                            WARN("Unknown error type at TRUTH %s:%d", ctg.data(), ctg_vars[TRUTH][h]->poss[i]);
                            continue;
                        }
                        // corresponding query call is only correct until its Qscore, after which it falls below
                        // the quality threshold, is filtered, and becomes a false negative
                        for (int qual = g.min_qual; qual <= q; qual++) {
                            truth_counts[t][ ctg_vars[TRUTH][h]->errtypes[swap][i] ][qual-g.min_qual]++;
                        }
                        for (int qual = q+1; qual < g.max_qual; qual++) {
                            truth_counts[t][ERRTYPE_FN][qual-g.min_qual]++;
                        }
                    }
                }
            }
        }
    }

    // write results
    std::string out_pr_fn = g.out_prefix + "precision-recall.tsv";
    FILE* out_pr = 0;
    if (g.write) {
        if (g.verbosity >= 1) INFO("  Writing precision-recall results to '%s'", out_pr_fn.data());
        out_pr = fopen(out_pr_fn.data(), "w");
        fprintf(out_pr, "VAR_TYPE\tMIN_QUAL\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE\t"
                "TRUTH_TOTAL\tTRUTH_TP\tTRUTH_FN\tQUERY_TOTAL\tQUERY_TP\tQUERY_FP\n");
    }
    std::vector<float> max_f1_score = {0, 0};
    std::vector<int> max_f1_qual = {0, 0};
    for (int type = 0; type < VARTYPES; type++) {

        // only sweeping query qualities; always consider all truth variants
        for (int qual = g.min_qual; qual <= g.max_qual; qual++) {
            int qidx = qual - g.min_qual;

            // define helper variables
            int query_tp = query_counts[type][ERRTYPE_TP][qidx];
            int query_fp = query_counts[type][ERRTYPE_FP][qidx];
            int query_tot = query_tp + query_fp;
            int truth_tp = truth_counts[type][ERRTYPE_TP][qidx];
            int truth_fn = truth_counts[type][ERRTYPE_FN][qidx];
            int truth_tot = truth_tp + truth_fn;

            // calculate summary metrics
            float precision = query_tot == 0 ? 1 : float(query_tp) / query_tot;
            float recall = truth_tot == 0 ? 1 : float(truth_tp) / truth_tot;
            float f1_score = precision+recall ? 2*precision*recall / (precision + recall) : 0;
            if (f1_score > max_f1_score[type]) {
                max_f1_score[type] = f1_score;
                max_f1_qual[type] = qual;
            }

            if (g.write) fprintf(out_pr, 
                    "%s\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    vartype_strs[type].data(),
                    qual,
                    precision,
                    recall, 
                    f1_score,
                    qscore(1-f1_score),
                    truth_tot,
                    int(truth_counts[type][ERRTYPE_TP][qidx]),
                    truth_fn,
                    query_tot, 
                    int(query_counts[type][ERRTYPE_TP][qidx]),
                    int(query_counts[type][ERRTYPE_FP][qidx])
           );
        }
    }
    if (g.write) fclose(out_pr);

    // print summary output
    std::string out_pr_summ_fn = g.out_prefix + "precision-recall-summary.tsv";
    FILE* out_pr_summ = 0;
    if (g.write) {
        if (g.verbosity >= 1) 
            INFO("  Writing precision-recall summary to '%s'", out_pr_summ_fn.data());
        out_pr_summ = fopen(out_pr_summ_fn.data(), "w");
        fprintf(out_pr_summ, "VAR_TYPE\tMIN_QUAL\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\t\tRECALL\t\tF1_SCORE\tF1_QSCORE\n");
    }
    INFO(" ");
    INFO("%sPRECISION-RECALL SUMMARY%s", COLOR_BLUE, COLOR_WHITE);
    for (int type = 0; type < VARTYPES; type++) {
        std::vector<int> quals = {g.min_qual, max_f1_qual[type]};
        INFO(" ");
        INFO("%sTYPE\tMIN_QUAL\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\t\tRECALL\t\tF1_SCORE\tF1_QSCORE%s",
                COLOR_BLUE, COLOR_WHITE);

        for (int qual : quals) {
            // redo calculations for these two
            int qidx = qual - g.min_qual;

            // define helper variables
            int query_tp = query_counts[type][ERRTYPE_TP][qidx];
            int query_fp = query_counts[type][ERRTYPE_FP][qidx];
            int query_tot = query_tp + query_fp;
            if (query_tot == 0) WARN("No QUERY %s variant calls.", vartype_strs[type].data());
            int truth_tp = truth_counts[type][ERRTYPE_TP][qidx];
            int truth_fn = truth_counts[type][ERRTYPE_FN][qidx];
            int truth_tot = truth_tp + truth_fn;
            if (truth_tot == 0) WARN("No TRUTH %s variant calls.", vartype_strs[type].data());

            // calculate summary metrics
            float precision = query_tot == 0 ? 1 : float(query_tp) / query_tot;
            float recall = truth_tot == 0 ? 1 : float(truth_tp) / truth_tot;
            float f1_score = precision+recall ? 2*precision*recall / (precision + recall) : 0;

            // print summary
            INFO("%s%s\tQ >= %d\t\t%-16d%-16d%-16d%-16d%f\t%f\t%f\t%f%s",
                COLOR_BLUE,
                vartype_strs[type].data(),
                qual,
                int(truth_counts[type][ERRTYPE_TP][qidx]),
                int(query_counts[type][ERRTYPE_TP][qidx]),
                truth_fn,
                int(query_counts[type][ERRTYPE_FP][qidx]),
                precision,
                recall, 
                f1_score,
                qscore(1-f1_score),
                COLOR_WHITE
            );
            if (g.write) fprintf(out_pr_summ,
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
    if (g.write) fclose(out_pr_summ);
}


/*******************************************************************************/


void write_distance(const editData & edits) {

    // log all distance results
    std::string dist_fn = g.out_prefix + "distance.tsv";
    FILE* out_dists = 0;
    if (g.write) {
        out_dists = fopen(dist_fn.data(), "w");
        fprintf(out_dists, "MIN_QUAL\tSUB_DE\tINS_DE\tDEL_DE\tSUB_ED\tINS_ED\tDEL_ED\t"
                "DISTINCT_EDITS\tEDIT_DIST\tALN_SCORE\tALN_QSCORE\n");
        if (g.verbosity >= 1) INFO("  Writing distance results to '%s'", dist_fn.data());
    }

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

        if (g.write) fprintf(out_dists, 
                "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", 
                q, distinct_edits[TYPE_SUB],
                distinct_edits[TYPE_INS], distinct_edits[TYPE_DEL],
                edit_dists[TYPE_SUB], edit_dists[TYPE_INS], edit_dists[TYPE_DEL],
                distinct_edits[TYPE_ALL], edit_dists[TYPE_ALL],
                edits.get_score(q), qscore(double(edits.get_score(q))/orig_score));
    }
    if (g.write) fclose(out_dists);

    // summarize distance results
    std::string dist_summ_fn = g.out_prefix + "distance-summary.tsv";
    FILE* dists_summ = 0;
    if (g.write) {
        dists_summ = fopen(dist_summ_fn.data(), "w");
        if (g.verbosity >= 1) 
            INFO("  Writing distance summary to '%s'", dist_summ_fn.data());
        fprintf(dists_summ, "VAR_TYPE\tMIN_QUAL\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE\tALN_QSCORE\n");
    }
    INFO(" ");
    INFO("%sALIGNMENT DISTANCE SUMMARY%s", COLOR_BLUE, COLOR_WHITE);
    for (int type = 0; type < TYPES; type++) {

        // skip INS/DEL individually unless higher print verbosity
        if (g.verbosity == 0 && type != TYPE_ALL)
            continue;
        if (g.verbosity == 1 && (type == TYPE_INS || type == TYPE_DEL))
            continue;

        INFO(" ");
        if (type == TYPE_ALL) {
            INFO("%sTYPE\tMIN_QUAL\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE\tALN_QSCORE%s",
                    COLOR_BLUE, COLOR_WHITE);
        } else {
            INFO("%sTYPE\tMIN_QUAL\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE%s",
                    COLOR_BLUE, COLOR_WHITE);
        }
        std::vector<int> quals = {g.min_qual, best_qual[type], g.max_qual+1};
        for (int q : quals) {

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
            if (g.write) fprintf(dists_summ, "%s\t%d\t%d\t%d\t%f\t%f\t%f\n", type_strs2[type].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, all_qscore);
            if (type == TYPE_ALL) {
                INFO("%s%s\tQ >= %d\t\t%-16d%-16d%f\t%f\t%f%s", COLOR_BLUE, type_strs2[type].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, all_qscore, COLOR_WHITE);
            } else {
                INFO("%s%s\tQ >= %d\t\t%-16d%-16d%f\t%f%s", COLOR_BLUE, type_strs2[type].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, COLOR_WHITE);
            }
        }
    }
    INFO(" ");
    if (g.write) fclose(dists_summ);
}


/*******************************************************************************/


void write_results(
        std::unique_ptr<phaseblockData> & phasedata_ptr, 
        const editData & edits) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[8/8] Writing results%s", COLOR_PURPLE, COLOR_WHITE);

    // print summary (precision/recall) information
    write_precision_recall(phasedata_ptr);

    // print distance information
    write_distance(edits);

    // print edit information
    if (g.write) {
        std::string edit_fn = g.out_prefix + "edits.tsv";
        FILE* out_edits = fopen(edit_fn.data(), "w");
        fprintf(out_edits, "CONTIG\tSTART\tHAP\tTYPE\tSIZE\tSUPERCLUSTER\tMIN_QUAL\tMAX_QUAL\n");
        if (g.verbosity >= 1) INFO("  Writing edit results to '%s'", edit_fn.data());
        for (int i = 0; i < edits.n; i++) {
            fprintf(out_edits, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", 
                    edits.ctgs[i].data(), edits.poss[i], edits.haps[i],
                    type_strs[edits.types[i]].data(), edits.lens[i],
                    edits.superclusters[i], edits.min_quals[i], edits.max_quals[i]);
        }
        fclose(out_edits);

        // print phasing information
        std::string out_phaseblocks_fn = g.out_prefix + "phase-blocks.tsv";
        FILE* out_phaseblocks = fopen(out_phaseblocks_fn.data(), "w");
        if (g.verbosity >= 1) INFO("  Writing phasing results to '%s'", out_phaseblocks_fn.data());
        fprintf(out_phaseblocks, "CONTIG\tPHASE_BLOCK\tSTART\tSTOP\tSIZE\tSUPERCLUSTERS\tBLOCK_STATE\n");
        for (std::string ctg : phasedata_ptr->contigs) {
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            for (int i = 0; i < ctg_pbs->n && ctg_scs->n > 0; i++) {
                int beg_idx = ctg_pbs->phase_blocks[i];
                int end_idx = ctg_pbs->phase_blocks[i+1]-1;
                int beg = ctg_scs->begs[beg_idx];
                int end = ctg_scs->ends[end_idx];
                int phase = ctg_pbs->block_states[i];
                fprintf(out_phaseblocks, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
                        ctg.data(), i, beg, end, end-beg, end_idx-beg_idx+1, phase);
            }
        }
        fclose(out_phaseblocks);

        // print clustering information
        std::string out_clusterings_fn = g.out_prefix + "superclusters.tsv";
        FILE* out_clusterings = fopen(out_clusterings_fn.data(), "w");
        if (g.verbosity >= 1) INFO("  Writing superclustering results to '%s'", out_clusterings_fn.data());
        fprintf(out_clusterings, "CONTIG\tSUPERCLUSTER\tSTART\tSTOP\tSIZE\tQUERY1_VARS\tQUERY2_VARS"
                "\tTRUTH1_VARS\tTRUTH2_VARS\tORIG_ED\tSWAP_ED\tPHASE\tPHASE_SET\tPHASE_BLOCK\tFLIP_ERROR\n");
        for (std::string ctg : phasedata_ptr->contigs) {
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            int phase_block_idx = 0;
            for (int i = 0; i < ctg_scs->n; i++) {

                // we've entered the next phase block
                if (i >= ctg_pbs->phase_blocks[phase_block_idx+1]) {
                    phase_block_idx++;
                }

                // set supercluster flip/swap based on phaseblock and sc phasing
                bool phase_switch = ctg_pbs->block_states[phase_block_idx];
                int phase_sc = ctg_scs->phase[i];
                bool flip_error;
                if (phase_switch) {
                    flip_error = phase_sc == PHASE_ORIG;
                } else { // no phase switch
                    flip_error = phase_sc == PHASE_SWAP;
                }

                // count query vars, allowing empty haps
                int query1_vars = ctg_scs->ctg_variants[QUERY][HAP1]->clusters.size() ?
                    ctg_scs->ctg_variants[QUERY][HAP1]->clusters[
                            ctg_scs->superclusters[QUERY][HAP1][i+1]] -
                        ctg_scs->ctg_variants[QUERY][HAP1]->clusters[
                            ctg_scs->superclusters[QUERY][HAP1][i]] : 0;
                int query2_vars = ctg_scs->ctg_variants[QUERY][HAP2]->clusters.size() ?
                    ctg_scs->ctg_variants[QUERY][HAP2]->clusters[
                            ctg_scs->superclusters[QUERY][HAP2][i+1]] -
                        ctg_scs->ctg_variants[QUERY][HAP2]->clusters[
                            ctg_scs->superclusters[QUERY][HAP2][i]] : 0;
                int truth1_vars = ctg_scs->ctg_variants[TRUTH][HAP1]->clusters.size() ?
                    ctg_scs->ctg_variants[TRUTH][HAP1]->clusters[
                            ctg_scs->superclusters[TRUTH][HAP1][i+1]] -
                        ctg_scs->ctg_variants[TRUTH][HAP1]->clusters[
                            ctg_scs->superclusters[TRUTH][HAP1][i]] : 0;
                int truth2_vars = ctg_scs->ctg_variants[TRUTH][HAP2]->clusters.size() ?
                    ctg_scs->ctg_variants[TRUTH][HAP2]->clusters[
                            ctg_scs->superclusters[TRUTH][HAP2][i+1]] -
                        ctg_scs->ctg_variants[TRUTH][HAP2]->clusters[
                            ctg_scs->superclusters[TRUTH][HAP2][i]] : 0;

                // print data
                fprintf(out_clusterings, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%s\n", 
                    ctg.data(), i, ctg_scs->begs[i], ctg_scs->ends[i],
                    ctg_scs->ends[i] - ctg_scs->begs[i],
                    query1_vars, query2_vars, truth1_vars, truth2_vars,
                    ctg_scs->orig_phase_dist[i], ctg_scs->swap_phase_dist[i],
                    phase_strs[phase_sc].data(), ctg_scs->phase_sets[i], phase_block_idx, flip_error ? "TRUE" : "FALSE"
               );
            }
        }
        fclose(out_clusterings);

        // print query variant information
        std::string out_query_fn = g.out_prefix + "query.tsv";
        if (g.verbosity >= 1) INFO("  Writing call variant results to '%s'", out_query_fn.data());
        FILE* out_query = fopen(out_query_fn.data(), "w");
        fprintf(out_query, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERR_TYPE"
                "\tCREDIT\tCLUSTER\tSUPERCLUSTER\tSYNC_GROUP\tLOCATION\n");
        for (std::string ctg : phasedata_ptr->contigs) {

            // set pointers to variants and superclusters
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            std::shared_ptr<ctgVariants> query1_vars = ctg_scs->ctg_variants[QUERY][HAP1];
            std::shared_ptr<ctgVariants> query2_vars = ctg_scs->ctg_variants[QUERY][HAP2];

            int supercluster_idx = 0;
            int var1_idx = 0;
            int var2_idx = 0;
            int cluster1_idx = 0;
            int cluster2_idx = 0;
            int phase_block = 0;
            while (var1_idx < query1_vars->n || var2_idx < query2_vars->n) {
                if (var2_idx >= query2_vars->n || // only query1 has remaining vars
                        (var1_idx < query1_vars->n && query1_vars->poss[var1_idx] < query2_vars->poss[var2_idx])) { // query1 var next

                    // we've entered the next supercluster
                    if (cluster1_idx+1 >= int(query1_vars->clusters.size())) 
                        ERROR("Out of bounds cluster during write_results(): query1")
                    if (query1_vars->clusters[cluster1_idx+1] <= var1_idx) cluster1_idx++;
                    if (supercluster_idx >= int(ctg_scs->begs.size())) 
                        ERROR("Out of bounds supercluster during write_results(): query1")
                    while (query1_vars->poss[var1_idx] >= ctg_scs->ends[supercluster_idx])
                        supercluster_idx++;

                    // update phasing info
                    if (supercluster_idx >= ctg_pbs->phase_blocks[phase_block+1])
                        phase_block++;
                    bool phase_switch = ctg_pbs->block_states[phase_block];
                    int phase_sc = ctg_scs->phase[supercluster_idx];
                    bool swap;
                    switch (phase_sc) {
                        case PHASE_ORIG: swap = false; break;
                        case PHASE_SWAP: swap = true; break;
                        default: swap = phase_switch; break;
                    }

                    fprintf(out_query, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%d\t%s\n",
                            ctg.data(),
                            query1_vars->poss[var1_idx],
                            query1_vars->haps[var1_idx],
                            query1_vars->refs[var1_idx].data(),
                            query1_vars->alts[var1_idx].data(),
                            query1_vars->var_quals[var1_idx],
                            type_strs[query1_vars->types[var1_idx]].data(),
                            error_strs[query1_vars->errtypes[swap][var1_idx]].data(),
                            query1_vars->credit[swap][var1_idx],
                            cluster1_idx,
                            supercluster_idx,
                            query1_vars->sync_group[swap][var1_idx],
                            region_strs[query1_vars->locs[var1_idx]].data()
                           );
                    var1_idx++;
                } else { // process query2 var
                    if (cluster2_idx+1 >= int(query2_vars->clusters.size()))
                        ERROR("Out of bounds cluster during write_results(): query2")
                    if (query2_vars->clusters[cluster2_idx+1] <= var2_idx) cluster2_idx++;
                    if (supercluster_idx >= int(ctg_scs->begs.size())) 
                        ERROR("Out of bounds supercluster during write_results(): query2")
                    while (query2_vars->poss[var2_idx] >= ctg_scs->ends[supercluster_idx])
                        supercluster_idx++;

                    // update phasing info
                    if (supercluster_idx >= ctg_pbs->phase_blocks[phase_block+1])
                        phase_block++;
                    bool phase_switch = ctg_pbs->block_states[phase_block];
                    int phase_sc = ctg_scs->phase[supercluster_idx];
                    bool swap;
                    switch (phase_sc) {
                        case PHASE_ORIG: swap = false; break;
                        case PHASE_SWAP: swap = true; break;
                        default: swap = phase_switch; break;
                    }

                    fprintf(out_query, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%d\t%s\n",
                            ctg.data(),
                            query2_vars->poss[var2_idx],
                            query2_vars->haps[var2_idx],
                            query2_vars->refs[var2_idx].data(),
                            query2_vars->alts[var2_idx].data(),
                            query2_vars->var_quals[var2_idx],
                            type_strs[query2_vars->types[var2_idx]].data(),
                            error_strs[query2_vars->errtypes[swap][var2_idx]].data(),
                            query2_vars->credit[swap][var2_idx],
                            cluster2_idx,
                            supercluster_idx,
                            query2_vars->sync_group[swap][var2_idx],
                            region_strs[query2_vars->locs[var2_idx]].data()
                           );
                    var2_idx++;
                }
            }
        }
        fclose(out_query);
        
        // print truth variant information
        std::string out_truth_fn = g.out_prefix + "truth.tsv";
        if (g.verbosity >= 1) INFO("  Writing truth variant results to '%s'", out_truth_fn.data());
        FILE* out_truth = fopen(out_truth_fn.data(), "w");
        fprintf(out_truth, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERRTYPE\tCREDIT\tCLUSTER\tSUPERCLUSTER\tSYNC_GROUP\tLOCATION\n");
        for (std::string ctg : phasedata_ptr->contigs) {

            // set pointers to variants and superclusters
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            std::shared_ptr<ctgVariants> truth1_vars = ctg_scs->ctg_variants[TRUTH][HAP1];
            std::shared_ptr<ctgVariants> truth2_vars = ctg_scs->ctg_variants[TRUTH][HAP2];

            int var1_idx = 0;
            int var2_idx = 0;
            int cluster1_idx = 0;
            int cluster2_idx = 0;
            int supercluster_idx = 0;
            int phase_block = 0;
            while (var1_idx < truth1_vars->n || var2_idx < truth2_vars->n) {
                if (var2_idx >= truth2_vars->n || // only truth1 has remaining vars
                        (var1_idx < truth1_vars->n && truth1_vars->poss[var1_idx] < truth2_vars->poss[var2_idx])) { // truth1 var next
                    if (cluster1_idx+1 >= int(truth1_vars->clusters.size()))
                        ERROR("Out of bounds cluster during write_results(): truth1")
                    if (truth1_vars->clusters[cluster1_idx+1] <= var1_idx) cluster1_idx++;
                    if (supercluster_idx >= int(ctg_scs->begs.size())) 
                        ERROR("Out of bounds supercluster during write_results(): truth1")
                    while (truth1_vars->poss[var1_idx] >= ctg_scs->ends[supercluster_idx])
                        supercluster_idx++;

                    // update phasing info
                    if (supercluster_idx >= ctg_pbs->phase_blocks[phase_block+1])
                        phase_block++;
                    bool phase_switch = ctg_pbs->block_states[phase_block];
                    int phase_sc = ctg_scs->phase[supercluster_idx];
                    bool swap;
                    switch (phase_sc) {
                        case PHASE_ORIG: swap = false; break;
                        case PHASE_SWAP: swap = true; break;
                        default: swap = phase_switch; break;
                    }

                    fprintf(out_truth, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%d\t%s\n",
                            ctg.data(),
                            truth1_vars->poss[var1_idx],
                            truth1_vars->haps[var1_idx],
                            truth1_vars->refs[var1_idx].data(),
                            truth1_vars->alts[var1_idx].data(),
                            truth1_vars->var_quals[var1_idx],
                            type_strs[truth1_vars->types[var1_idx]].data(),
                            error_strs[truth1_vars->errtypes[swap][var1_idx]].data(),
                            truth1_vars->credit[swap][var1_idx],
                            cluster1_idx,
                            supercluster_idx,
                            truth1_vars->sync_group[swap][var1_idx],
                            region_strs[truth1_vars->locs[var1_idx]].data()
                           );
                    var1_idx++;
                } else { // process truth2 var
                    if (cluster2_idx+1 >= int(truth2_vars->clusters.size())) {
                        ERROR("Out of bounds cluster during write_results(): truth2")
                    }
                    if (truth2_vars->clusters[cluster2_idx+1] <= var2_idx) cluster2_idx++;
                    if (supercluster_idx >= int(ctg_scs->begs.size())) 
                        ERROR("Out of bounds supercluster during write_results(): truth2")
                    while (truth2_vars->poss[var2_idx] >= ctg_scs->ends[supercluster_idx])
                        supercluster_idx++;

                    // update phasing info
                    if (supercluster_idx >= ctg_pbs->phase_blocks[phase_block+1])
                        phase_block++;
                    bool phase_switch = ctg_pbs->block_states[phase_block];
                    int phase_sc = ctg_scs->phase[supercluster_idx];
                    bool swap;
                    switch (phase_sc) {
                        case PHASE_ORIG: swap = false; break;
                        case PHASE_SWAP: swap = true; break;
                        default: swap = phase_switch; break;
                    }

                    fprintf(out_truth, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%d\t%s\n",
                            ctg.data(),
                            truth2_vars->poss[var2_idx],
                            truth2_vars->haps[var2_idx],
                            truth2_vars->refs[var2_idx].data(),
                            truth2_vars->alts[var2_idx].data(),
                            truth2_vars->var_quals[var2_idx],
                            type_strs[truth2_vars->types[var2_idx]].data(),
                            error_strs[truth2_vars->errtypes[swap][var2_idx]].data(),
                            truth2_vars->credit[swap][var2_idx],
                            cluster2_idx,
                            supercluster_idx,
                            truth2_vars->sync_group[swap][var2_idx],
                            region_strs[truth2_vars->locs[var2_idx]].data()
                           );
                    var2_idx++;
                }
            }
        }
        fclose(out_truth);
    }
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
