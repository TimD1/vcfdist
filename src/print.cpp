#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

#include "print.h"
#include "dist.h"

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

float qscore(double p_error) {
    return std::min(100.0, std::max(0.0, -10 * std::log10(p_error)));
}

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
        "min_var_qual = %d\nmax_var_qual = %d\nmax_var_size = %d\nsv_threshold = %d\n"
        "phase_threshold = %f\ncredit_threshold = %f\nrealign_truth = %s\nrealign_query = %s\n"
        "realign_only = %s\ncluster_method = '%s'\ncluster_min_gap = %d\n"
        "reach_min_gap = %d\nmax_cluster_itrs = %d\nmax_threads = %d\nmax_ram = %f\n"
        "sub = %d\nopen = %d\nextend = %d\n",
        g.PROGRAM.data(), g.VERSION.data(), g.out_prefix.data(), g.cmd.data(), g.ref_fasta_fn.data(), 
        g.query_vcf_fn.data(), g.truth_vcf_fn.data(), g.bed_fn.data(), b2s(g.write).data(), 
        filters_str.data(), g.min_qual, g.max_qual, g.max_size, g.sv_threshold,
        g.phase_threshold, g.credit_threshold, b2s(g.realign_truth).data(), b2s(g.realign_query).data(),
        b2s(g.realign_only).data(), g.cluster_method.data(), g.cluster_min_gap,
        g.reach_min_gap, g.max_cluster_itrs, g.max_threads, g.max_ram,
        g.sub, g.open, g.extend);
    fclose(out_params);
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


void write_precision_recall(std::unique_ptr<phaseblockData> & phasedata_ptr) {

    // for each class, store variant counts above each quality threshold
    // init counters; ax0: SNP/INDEL/SV/ALL, ax1: TP,FP,FN ax2: QUAL
    std::vector< std::vector< std::vector<float> > > query_counts(VARTYPES,
            std::vector< std::vector<float> >(ERRTYPES, 
            std::vector<float>(g.max_qual-g.min_qual+1, 0.0))) ;
    std::vector< std::vector< std::vector<float> > > truth_counts(VARTYPES,
            std::vector< std::vector<float> >(ERRTYPES, 
            std::vector<float>(g.max_qual-g.min_qual+1, 0.0))) ;

    // calculate summary statistics
    for (std::string ctg : phasedata_ptr->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
        std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
        std::shared_ptr<ctgVariants> qvars = ctg_scs->callset_vars[QUERY];
        std::shared_ptr<ctgVariants> tvars = ctg_scs->callset_vars[TRUTH];

        // add query
        if (ctg_scs->n)
        for (int vi = 0; vi < qvars->n; vi++) {
            for (int hi = 0; hi < HAPS; hi++) {
                if (!qvars->var_on_hap(vi, hi, true)) continue;
                float q = qvars->callq[hi][vi];
                int t = 0;
                if (qvars->types[vi] == TYPE_SUB) { // SNP
                    t = VARTYPE_SNP;
                } else if ((qvars->types[vi] == TYPE_INS && // small INDEL
                            int(qvars->alts[vi].size()) < g.sv_threshold) ||
                        (qvars->types[vi] == TYPE_DEL &&
                         int(qvars->refs[vi].size()) < g.sv_threshold)) {
                    t = VARTYPE_INDEL;
                } else { // SV
                    t = VARTYPE_SV;
                }
                if (qvars->errtypes[hi][vi] == ERRTYPE_UN) {
                    WARN("Unknown error type at QUERY %s:%d", ctg.data(), qvars->poss[vi]);
                    continue;
                }
                for (int qual = g.min_qual; qual <= q; qual++) {
                    query_counts[t][ qvars->errtypes[hi][vi] ][qual-g.min_qual]++;
                    query_counts[VARTYPE_ALL][ qvars->errtypes[hi][vi] ][qual-g.min_qual]++;
                }
            }
        }

        // add truth
        if (ctg_scs->n)
        for (int vi = 0; vi < tvars->n; vi++) {
            for (int hi = 0; hi < HAPS; hi++) {
                if (!tvars->var_on_hap(vi, hi, false)) continue;
                float q = tvars->callq[hi][vi];
                int t = 0;
                if (tvars->types[vi] == TYPE_SUB) {
                    t = VARTYPE_SNP;
                } else if ((tvars->types[vi] == TYPE_INS &&
                            int(tvars->alts[vi].size()) < g.sv_threshold) ||
                        (tvars->types[vi] == TYPE_DEL && 
                         int(tvars->refs[vi].size()) < g.sv_threshold)) {
                    t = VARTYPE_INDEL;
                } else {
                    t = VARTYPE_SV;
                }
                if (tvars->errtypes[hi][vi] == ERRTYPE_UN) {
                    WARN("Unknown error type at TRUTH %s:%d", ctg.data(), tvars->poss[vi]);
                    continue;
                }
                // corresponding query call is only correct until its Qscore, after which it falls below
                // the quality threshold, is filtered, and becomes a false negative
                for (int qual = g.min_qual; qual <= q; qual++) {
                    truth_counts[t][ tvars->errtypes[hi][vi] ][qual-g.min_qual]++;
                    truth_counts[VARTYPE_ALL][ tvars->errtypes[hi][vi] ][qual-g.min_qual]++;
                }
                for (int qual = q+1; qual <= g.max_qual; qual++) {
                    truth_counts[t][ERRTYPE_FN][qual-g.min_qual]++;
                    truth_counts[VARTYPE_ALL][ERRTYPE_FN][qual-g.min_qual]++;
                }
            }
        }
    }

    // write results
    std::string out_pr_fn = g.out_prefix + "precision-recall.tsv";
    FILE* out_pr = 0;
    if (g.write) {
        if (g.verbosity >= 1) INFO(" ");
        if (g.verbosity >= 1) INFO("  Writing precision-recall results to '%s'", out_pr_fn.data());
        out_pr = fopen(out_pr_fn.data(), "w");
        fprintf(out_pr, "VAR_TYPE\tMIN_QUAL\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE\t"
                "TRUTH_TOTAL\tTRUTH_TP\tTRUTH_FN\tQUERY_TOTAL\tQUERY_TP\tQUERY_FP\n");
    }
    std::vector<float> max_f1_score(VARTYPES, 0);
    std::vector<int> max_f1_qual(VARTYPES, 0);
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
        fprintf(out_pr_summ, "VAR_TYPE\tTHRESHOLD\tMIN_QUAL\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\tRECALL\tF1_SCORE\tF1_QSCORE\n");
    }
    INFO(" ");
    INFO("%sPRECISION-RECALL SUMMARY%s", COLOR_BLUE, COLOR_WHITE);
    INFO(" ");
    INFO("%sTYPE\tTHRESHOLD\tTRUTH_TP\tQUERY_TP\tTRUTH_FN\tQUERY_FP\tPREC\t\tRECALL\t\tF1_SCORE\tF1_QSCORE%s",
            COLOR_BLUE, COLOR_WHITE);
    for (int type = 0; type < VARTYPES; type++) {
        std::vector<int> quals = {g.min_qual, max_f1_qual[type]};
        std::vector<std::string> thresholds = {"NONE", "BEST"};

        for (int i = 0; i < int(quals.size()); i++) {
            // redo calculations for these two
            int qual = quals[i];
            std::string thresh = thresholds[i];
            int qidx = qual - g.min_qual;

            // define helper variables
            int query_tp = query_counts[type][ERRTYPE_TP][qidx];
            int query_fp = query_counts[type][ERRTYPE_FP][qidx];
            int query_tot = query_tp + query_fp;
            if (query_tot == 0 && type != VARTYPE_ALL) 
                WARN("No QUERY %s variants pass all filters.", vartype_strs[type].data());
            int truth_tp = truth_counts[type][ERRTYPE_TP][qidx];
            int truth_fn = truth_counts[type][ERRTYPE_FN][qidx];
            int truth_tot = truth_tp + truth_fn;
            if (truth_tot == 0 && type != VARTYPE_ALL) 
                WARN("No TRUTH %s variants pass all filters.", vartype_strs[type].data());

            // calculate summary metrics
            float precision = query_tot == 0 ? 1 : float(query_tp) / query_tot;
            float recall = truth_tot == 0 ? 1 : float(truth_tp) / truth_tot;
            float f1_score = precision+recall > 0 ? 2*precision*recall / (precision + recall) : 0;

            // print summary
            INFO("%s%s\t%s Q >= %-2d\t%-16d%-16d%-16d%-16d%f\t%f\t%f\t%f%s",
                COLOR_BLUE,
                vartype_strs[type].data(),
                thresh.data(),
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
               "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",
                vartype_strs[type].data(),
                thresh.data(),
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
        INFO(" ");
    }
    if (g.write) fclose(out_pr_summ);
}


/*******************************************************************************/


void write_results(std::unique_ptr<phaseblockData> & phasedata_ptr) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[7/7] Writing results%s", COLOR_PURPLE, COLOR_WHITE);

    // print summary (precision/recall) information
    write_precision_recall(phasedata_ptr);

    if (g.write) {

        // print phasing information
        std::string out_phaseblocks_fn = g.out_prefix + "phase-blocks.tsv";
        FILE* out_phaseblocks = fopen(out_phaseblocks_fn.data(), "w");
        if (g.verbosity >= 1) INFO("  Writing phasing results to '%s'", out_phaseblocks_fn.data());
        fprintf(out_phaseblocks, "CONTIG\tPHASE_BLOCK\tSTART\tSTOP\tSIZE\tSUPERCLUSTERS\tFLIP_ERRORS\tSWITCH_ERRORS\n");
        for (std::string ctg : phasedata_ptr->contigs) {
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            for (int i = 0; i < ctg_pbs->n && ctg_scs->n > 0; i++) {
                int beg_idx = ctg_pbs->phase_blocks[i];
                int end_idx = ctg_pbs->phase_blocks[i+1]-1;
                int beg = ctg_scs->begs[beg_idx];
                int end = ctg_scs->ends[end_idx];
                int nswitches = 0;
                for (int si = 0; si < ctg_pbs->nswitches; si++) {
                    if (ctg_pbs->switches[si] > beg_idx && ctg_pbs->switches[si] <= end_idx) nswitches++;
                }
                int nflips = 0;
                for (int fi = 0; fi < ctg_pbs->nflips; fi++) {
                    if (ctg_pbs->flips[fi] > beg_idx && ctg_pbs->flips[fi] <= end_idx) nflips++;
                }
                fprintf(out_phaseblocks, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
                        ctg.data(), i, beg, end, end-beg, end_idx-beg_idx+1, nflips, nswitches);
            }
        }
        fclose(out_phaseblocks);

        // print clustering information
        std::string out_clusterings_fn = g.out_prefix + "superclusters.tsv";
        FILE* out_clusterings = fopen(out_clusterings_fn.data(), "w");
        if (g.verbosity >= 1) INFO("  Writing superclustering results to '%s'", out_clusterings_fn.data());
        fprintf(out_clusterings, "CONTIG\tSUPERCLUSTER\tSTART\tSTOP\tSIZE\t"
                "QUERY_VARS\tTRUTH_VARS\n");
        for (std::string ctg : phasedata_ptr->contigs) {
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            for (int i = 0; i < ctg_scs->n; i++) {

                // count query vars, allowing empty haps
                int query_vars = ctg_scs->callset_vars[QUERY]->clusters.size() ?
                    ctg_scs->callset_vars[QUERY]->clusters[
                            ctg_scs->superclusters[QUERY][i+1]] -
                        ctg_scs->callset_vars[QUERY]->clusters[
                            ctg_scs->superclusters[QUERY][i]] : 0;
                int truth_vars = ctg_scs->callset_vars[TRUTH]->clusters.size() ?
                    ctg_scs->callset_vars[TRUTH]->clusters[
                            ctg_scs->superclusters[TRUTH][i+1]] -
                        ctg_scs->callset_vars[TRUTH]->clusters[
                            ctg_scs->superclusters[TRUTH][i]] : 0;

                // print data
                fprintf(out_clusterings, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
                    ctg.data(), i, ctg_scs->begs[i], ctg_scs->ends[i],
                    ctg_scs->ends[i] - ctg_scs->begs[i], query_vars, truth_vars);
            }
        }
        fclose(out_clusterings);

        // print query variant information
        std::string out_query_fn = g.out_prefix + "query.tsv";
        if (g.verbosity >= 1) INFO("  Writing query variant results to '%s'", out_query_fn.data());
        FILE* out_query = fopen(out_query_fn.data(), "w");
        fprintf(out_query, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERR_TYPE"
                "\tCREDIT\tCLUSTER\tSUPERCLUSTER\tSYNC_GROUP\tREF_DIST\tQUERY_DIST\tLOCATION\n");
        for (std::string ctg : phasedata_ptr->contigs) {

            // set pointers to variants and superclusters
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            std::shared_ptr<ctgVariants> qvars = ctg_scs->callset_vars[QUERY];

            int sci = 0;
            int clust_idx = 0;
            for (int vi = 0; vi < qvars->n; vi++) {
                for (int hi = 0; hi < HAPS; hi++) {
                    if (!qvars->var_on_hap(vi, hi, /*calc=*/ true)) continue;

                    // update cluster and supercluster
                    if (clust_idx+1 >= int(qvars->clusters.size())) 
                        ERROR("Out of bounds cluster during write_results(): query")
                    if (qvars->clusters[clust_idx+1] <= vi) clust_idx++;
                    if (sci >= int(ctg_scs->begs.size())) 
                        ERROR("Out of bounds supercluster during write_results(): query")
                    while (qvars->poss[vi] >= ctg_scs->ends[sci])
                        sci++;

                    fprintf(out_query, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%s\n",
                            ctg.data(),
                            qvars->poss[vi],
                            hi,
                            qvars->refs[vi].data(),
                            qvars->alts[vi].data(),
                            qvars->var_quals[vi],
                            type_strs[qvars->types[vi]].data(),
                            error_strs[qvars->errtypes[hi][vi]].data(),
                            qvars->credit[hi][vi],
                            clust_idx,
                            sci,
                            qvars->sync_group[hi][vi],
                            qvars->ref_ed[hi][vi],
                            qvars->query_ed[hi][vi],
                            region_strs[qvars->locs[vi]].data()
                           );
                }
            }
        }
        fclose(out_query);
        
        // print truth variant information
        std::string out_truth_fn = g.out_prefix + "truth.tsv";
        if (g.verbosity >= 1) INFO("  Writing truth variant results to '%s'", out_truth_fn.data());
        FILE* out_truth = fopen(out_truth_fn.data(), "w");
        fprintf(out_truth, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERRTYPE\tCREDIT\tCLUSTER\tSUPERCLUSTER\tSYNC_GROUP\tREF_DIST\tQUERY_DIST\tLOCATION\n");
        for (std::string ctg : phasedata_ptr->contigs) {

            // set pointers to variants and superclusters
            std::shared_ptr<ctgPhaseblocks> ctg_pbs = phasedata_ptr->phase_blocks[ctg];
            std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
            std::shared_ptr<ctgVariants> tvars = ctg_scs->callset_vars[TRUTH];

            int sci = 0;
            int clust_idx = 0;
            for (int vi = 0; vi < tvars->n; vi++) {
                for (int hi = 0; hi < HAPS; hi++) {
                    if (!tvars->var_on_hap(vi, hi, /*calc=*/ false)) continue;

                    // update cluster and supercluster
                    if (clust_idx+1 >= int(tvars->clusters.size()))
                        ERROR("Out of bounds cluster during write_results(): truth")
                    if (tvars->clusters[clust_idx+1] <= vi) clust_idx++;
                    if (sci >= int(ctg_scs->begs.size())) 
                        ERROR("Out of bounds supercluster during write_results(): truth")
                    while (tvars->poss[vi] >= ctg_scs->ends[sci])
                        sci++;

                    fprintf(out_truth, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%s\n",
                            ctg.data(),
                            tvars->poss[vi],
                            hi,
                            tvars->refs[vi].data(),
                            tvars->alts[vi].data(),
                            tvars->var_quals[vi],
                            type_strs[tvars->types[vi]].data(),
                            error_strs[tvars->errtypes[hi][vi]].data(),
                            tvars->credit[hi][vi],
                            clust_idx,
                            sci,
                            tvars->sync_group[hi][vi],
                            tvars->ref_ed[hi][vi],
                            tvars->query_ed[hi][vi],
                            region_strs[tvars->locs[vi]].data()
                           );
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


/*******************************************************************************/

/* Helper function for printing a query graph alignment.
 */ 
std::string get_ptr_repr(idx4 cell, const std::unordered_map<idx4,idx4> & ptrs) {
    if (ptrs.find(cell) == ptrs.end()) return "  .";
    idx4 prev = ptrs.at(cell);
    if (cell.qni == prev.qni && cell.tni == prev.tni) { // same matrix
        if (cell.qi == prev.qi+1 && cell.ti == prev.ti) { return "  |"; }
        else if (cell.qi == prev.qi && cell.ti == prev.ti+1) { return "  _"; }
        else if (cell.qi == prev.qi+1 && cell.ti == prev.ti+1) { return "  \\"; }
        else { return "?1"; } // invalid
    } else if (cell.qi == 0 && prev.qni < cell.qni) { // different query node, print node id
        std::ostringstream ostr;
        ostr << std::setfill('0') << std::setw(2) << (prev.qni);
        return "^" + ostr.str();
    } else if (cell.ti == 0 && prev.tni < cell.tni) { // different truth node, print node id
        std::ostringstream ostr;
        ostr << std::setfill('0') << std::setw(2) << (prev.tni);
        return "<" + ostr.str();
    } else { return "?2"; } // invalid
}

void print_graph_ptrs(const std::shared_ptr<Graph> graph,
        const std::unordered_map<idx4,idx4> & ptrs) {
    printf("    TRUTH");

    // print truth node indices
    int tni = 0;
    while (tni < graph->tnodes) {
        for (int ti = 0; ti < int(graph->tseqs[tni].length()); ti++) {
            if (ti == 0) {
                    std::ostringstream ostr;
                    ostr << std::setfill(' ') << std::setw(2) << tni;
                    printf("%s %s", ostr.str().data(), type_strs[graph->ttypes[tni]].data());
            } else {
                printf("   ");
            }
        }
        tni++;
    }
    printf("\n");

    // print truth graph sequence
    printf(" QUERY    ");
    tni = 0;
    while (tni < graph->tnodes) {
        for (int ti = 0; ti < int(graph->tseqs[tni].length()); ti++) {
            printf("  %c", graph->tseqs[tni][ti]);   
        }
        printf("   ");
        tni++;
    }
    printf("\n");

    // print matrix
    int qni = 0;
    int qi = 0;
    while (qni < graph->qnodes) {
        tni = -1;
        while (tni < graph->tnodes) {
            if (tni < 0) { // print query seqeunce
                std::ostringstream ostr;
                ostr << std::setfill(' ') << std::setw(2) << qni;
                printf("%s %s %c  ", 
                        ostr.str().data(),
                        type_strs[graph->qtypes[qni]].data(),
                        graph->qseqs[qni][qi]);
            } else { // print cell
                for (int ti = 0; ti < int(graph->tseqs[tni].length()); ti++) {
                    idx4 cell(qni,tni,qi,ti);
                    std::string s = get_ptr_repr(cell, ptrs);
                    printf("%s", s.data());
                }
                printf("   ");
            }
            tni++;
        }
        printf("\n");
        qi++;
        if (qi == int(graph->qseqs[qni].size())) { // start next query node
            qi = 0;
            qni++;
            printf("\n");
        }
    }
}
