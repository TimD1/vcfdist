#include <string>
#include <vector>

#include "print.h"
#include "dist.h"

std::string GREEN(char c) { return "\033[32m" + std::string(1,c) + "\033[0m"; }
std::string GREEN(std::string str) { return "\033[32m" + str + "\033[0m"; }
std::string RED(char c) { return "\033[31m" + std::string(1,c) + "\033[0m"; }
std::string RED(std::string str) { return "\033[31m" + str + "\033[0m"; }
std::string BLUE(int i) { return "\033[34m" + std::to_string(i) + "\033[0m"; }
std::string BLUE(char c) { return "\033[34m" + std::string(1,c) + "\033[0m"; }
std::string BLUE(std::string str) { return "\033[34m" + str + "\033[0m"; }
std::string YELLOW(char c) { return "\033[33m" + std::string(1,c) + "\033[0m"; }
std::string YELLOW(std::string str) { return "\033[33m" + str + "\033[0m"; }
std::string PURPLE(char c) { return "\033[35m" + std::string(1,c) + "\033[0m"; }
std::string PURPLE(std::string str) { return "\033[35m" + str + "\033[0m"; }

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
            else if (ptrs[alt_idx][ref_idx] & PTR_SWAP)
                ptr_str[alt_idx*2+1][ref_idx*2+1] = '+';
            else
                ptr_str[alt_idx*2+1][ref_idx*2+1] = '*';

            // movements
            if (ptrs[alt_idx][ref_idx] & PTR_DIAG)
                ptr_str[alt_idx*2][ref_idx*2] = '\\';
            if (ptrs[alt_idx][ref_idx] & PTR_SUB)
                ptr_str[alt_idx*2][ref_idx*2] = 'X';
            if (ptrs[alt_idx][ref_idx] & PTR_LEFT)
                ptr_str[alt_idx*2+1][ref_idx*2] = '-';
            if (ptrs[alt_idx][ref_idx] & PTR_UP)
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



void write_results(std::unique_ptr<phaseData> & phasedata_ptr) {
    INFO(" ");
    INFO("7. Writing results");

    // print phasing information
    std::string out_phasings_fn = g.out_prefix + "phaseblocks.tsv";
    FILE* out_phasings = fopen(out_phasings_fn.data(), "w");
    INFO("  Printing phasing results to '%s'", out_phasings_fn.data());
    fprintf(out_phasings, "CONTIG\tSTART\tSTOP\tSIZE\n");
    for (auto ctg_name : phasedata_ptr->contigs) {
        ctgPhasings ctg_phasings = phasedata_ptr->ctg_phasings[ctg_name];
        std::shared_ptr<ctgClusters> ctg_superclusters = ctg_phasings.ctg_superclusters;
        for (int i = 0; i <= ctg_phasings.nswitches && ctg_phasings.nswitches > 0; i++) {
            int beg_idx = ctg_phasings.phase_blocks[i];
            int end_idx = ctg_phasings.phase_blocks[i+1]-1;
            int beg = ctg_superclusters->begs[beg_idx];
            int end = ctg_superclusters->ends[end_idx];
            fprintf(out_phasings, "%s\t%d\t%d\t%d\n", 
                    ctg_name.data(), beg, end, end-beg);
        }
    }
    fclose(out_phasings);

    // print clustering information
    std::string out_clusterings_fn = g.out_prefix + "superclusters.tsv";
    FILE* out_clusterings = fopen(out_clusterings_fn.data(), "w");
    INFO("  Printing superclustering results to '%s'", out_clusterings_fn.data());
    fprintf(out_clusterings, "CONTIG\tSTART\tSTOP\tSIZE\tCALLS1_VARS\tCALLS2_VARS\tTRUTH1_VARS\tTRUTH2_VARS\tORIG_ED\tSWAP_ED\tPHASE\tPHASE_BLOCK\n");
    for (auto ctg_name : phasedata_ptr->contigs) {
        auto ctg_phasings = phasedata_ptr->ctg_phasings[ctg_name];
        std::shared_ptr<ctgClusters> ctg_supclusts = ctg_phasings.ctg_superclusters;
        int phase_block_idx = 0;
        for (int i = 0; i < ctg_supclusts->n; i++) {

            // we've entered the next phase block
            if (i >= ctg_phasings.phase_blocks[phase_block_idx]) {
                phase_block_idx++;
            }

            // print data
            fprintf(out_phasings, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\n", 
                ctg_name.data(), 
                ctg_supclusts->begs[i],
                ctg_supclusts->ends[i],
                ctg_supclusts->ends[i] - ctg_supclusts->begs[i],
                ctg_supclusts->calls1_end_idx[i] - ctg_supclusts->calls1_beg_idx[i],
                ctg_supclusts->calls2_end_idx[i] - ctg_supclusts->calls2_beg_idx[i],
                ctg_supclusts->truth1_end_idx[i] - ctg_supclusts->truth1_beg_idx[i],
                ctg_supclusts->truth2_end_idx[i] - ctg_supclusts->truth2_beg_idx[i],
                ctg_supclusts->orig_phase_dist[i],
                ctg_supclusts->swap_phase_dist[i],
                phase_strs[ctg_supclusts->phase[i]].data(),
                phase_block_idx
           );
        }
    }
    fclose(out_clusterings);

    // print calls clusters information
    std::string out_calls_clust_fn = g.out_prefix + "calls-clusters.tsv";
    INFO("  Printing calls clustering results to '%s'", out_calls_clust_fn.data());
    FILE* out_calls_clust = fopen(out_calls_clust_fn.data(), "w");
    fprintf(out_calls_clust, "CONTIG\tHAP\tSTART\tSTOP\tSIZE\tVARIANTS\tSUPERCLUSTER\n");
    for (auto ctg_name : phasedata_ptr->contigs) {

        // set pointers to variants and superclusters
        std::shared_ptr<ctgVariants> calls1_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->calls1_vars;
        std::shared_ptr<ctgVariants> calls2_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->calls2_vars;
        std::shared_ptr<ctgClusters> ctg_supclusts = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters;

        int supercluster_idx = 0;
        for(int clust_idx = 0; clust_idx < int(calls1_vars->clusters.size())-1; clust_idx++) {
            int beg_idx = calls1_vars->clusters[clust_idx];
            int end_idx = calls1_vars->clusters[clust_idx+1];
            int beg = calls1_vars->poss[beg_idx]-1;
            int end = calls1_vars->poss[end_idx-1] + calls1_vars->rlens[end_idx-1]+1;
            fprintf(out_calls_clust, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    ctg_name.data(), 0, beg, end, end-beg,
                    end_idx-beg_idx, supercluster_idx);
            if (beg >= ctg_supclusts->begs[supercluster_idx]) {
                supercluster_idx++;
            }
        }
        supercluster_idx = 0;
        for(int clust_idx = 0; clust_idx < int(calls2_vars->clusters.size())-1; clust_idx++) {
            int beg_idx = calls2_vars->clusters[clust_idx];
            int end_idx = calls2_vars->clusters[clust_idx+1];
            int beg = calls2_vars->poss[beg_idx]-1;
            int end = calls2_vars->poss[end_idx-1] + calls2_vars->rlens[end_idx-1]+1;
            fprintf(out_calls_clust, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    ctg_name.data(), 1, beg, end, end-beg,
                    end_idx-beg_idx, supercluster_idx);
            if (beg >= ctg_supclusts->begs[supercluster_idx]) {
                supercluster_idx++;
            }
        }
    }
    fclose(out_calls_clust);

    // print truth clusters information
    std::string out_truth_clust_fn = g.out_prefix + "truth-clusters.tsv";
    INFO("  Printing truth clustering results to '%s'", out_truth_clust_fn.data());
    FILE* out_truth_clust = fopen(out_truth_clust_fn.data(), "w");
    fprintf(out_truth_clust, "CONTIG\tHAP\tSTART\tSTOP\tSIZE\tVARIANTS\tSUPERCLUSTER\n");
    for (auto ctg_name : phasedata_ptr->contigs) {

        // set pointers to variants and superclusters
        std::shared_ptr<ctgVariants> truth1_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->truth1_vars;
        std::shared_ptr<ctgVariants> truth2_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->truth2_vars;
        std::shared_ptr<ctgClusters> ctg_supclusts = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters;

        int supercluster_idx = 0;
        for(int clust_idx = 0; clust_idx < int(truth1_vars->clusters.size())-1; clust_idx++) {
            int beg_idx = truth1_vars->clusters[clust_idx];
            int end_idx = truth1_vars->clusters[clust_idx+1];
            int beg = truth1_vars->poss[beg_idx]-1;
            int end = truth1_vars->poss[end_idx-1] + truth1_vars->rlens[end_idx-1]+1;
            fprintf(out_truth_clust, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    ctg_name.data(), 0, beg, end, end-beg,
                    end_idx-beg_idx, supercluster_idx);
            if (beg >= ctg_supclusts->begs[supercluster_idx]) {
                supercluster_idx++;
            }
        }
        supercluster_idx = 0;
        for(int clust_idx = 0; clust_idx < int(truth2_vars->clusters.size())-1; clust_idx++) {
            int beg_idx = truth2_vars->clusters[clust_idx];
            int end_idx = truth2_vars->clusters[clust_idx+1];
            int beg = truth2_vars->poss[beg_idx]-1;
            int end = truth2_vars->poss[end_idx-1] + truth2_vars->rlens[end_idx-1]+1;
            fprintf(out_truth_clust, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    ctg_name.data(), 1, beg, end, end-beg,
                    end_idx-beg_idx, supercluster_idx);
            if (beg >= ctg_supclusts->begs[supercluster_idx]) {
                supercluster_idx++;
            }
        }
    }
    fclose(out_truth_clust);

    // print calls variant information
    std::string out_calls_fn = g.out_prefix + "calls.tsv";
    INFO("  Printing call variant results to '%s'", out_calls_fn.data());
    FILE* out_calls = fopen(out_calls_fn.data(), "w");
    fprintf(out_calls, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERRTYPE\tCREDIT\tORIG_GT\tCLUSTER\tSUPERCLUSTER\tLOCATION\n");
    for (auto ctg_name : phasedata_ptr->contigs) {

        // set pointers to variants and superclusters
        std::shared_ptr<ctgVariants> calls1_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->calls1_vars;
        std::shared_ptr<ctgVariants> calls2_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->calls2_vars;
        std::shared_ptr<ctgClusters> ctg_supclusts = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters;

        int supercluster_idx = 0;
        int var1_idx = 0;
        int var2_idx = 0;
        int cluster1_idx = 0;
        int cluster2_idx = 0;
        while (var1_idx < calls1_vars->n || var2_idx < calls2_vars->n) {
            if (var2_idx >= calls2_vars->n || // only calls1 has remaining vars
                    calls1_vars->poss[var1_idx] < calls2_vars->poss[var2_idx]) { // calls1 var next

                // we've entered the next supercluster
                if (calls1_vars->clusters[cluster1_idx] <= var1_idx) cluster1_idx++;
                fprintf(out_calls, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%s\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        calls1_vars->poss[var1_idx],
                        calls1_vars->haps[var1_idx],
                        calls1_vars->refs[var1_idx].data(),
                        calls1_vars->alts[var1_idx].data(),
                        calls1_vars->var_quals[var1_idx],
                        type_strs[calls1_vars->types[var1_idx]].data(),
                        error_strs[calls1_vars->errtypes[var1_idx]].data(),
                        calls1_vars->credit[var1_idx],
                        gt_strs[calls1_vars->orig_gts[var1_idx]].data(),
                        cluster1_idx,
                        supercluster_idx,
                        region_strs[calls1_vars->locs[var1_idx]].data()
                       );
                if (calls1_vars->poss[var1_idx] >= ctg_supclusts->begs[supercluster_idx]) {
                    supercluster_idx++;
                }
                if (cluster1_idx >= int(calls1_vars->clusters.size())) 
                    ERROR("Out of bounds cluster during write_results(): calls1")
                var1_idx++;
            } else { // process calls2 var
                if (calls2_vars->clusters[cluster2_idx] <= var2_idx) cluster2_idx++;
                fprintf(out_calls, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%s\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        calls2_vars->poss[var2_idx],
                        calls2_vars->haps[var2_idx],
                        calls2_vars->refs[var2_idx].data(),
                        calls2_vars->alts[var2_idx].data(),
                        calls2_vars->var_quals[var2_idx],
                        type_strs[calls2_vars->types[var2_idx]].data(),
                        error_strs[calls2_vars->errtypes[var2_idx]].data(),
                        calls2_vars->credit[var2_idx],
                        gt_strs[calls2_vars->orig_gts[var2_idx]].data(),
                        cluster2_idx,
                        supercluster_idx,
                        region_strs[calls2_vars->locs[var2_idx]].data()
                       );
                if (calls2_vars->poss[var2_idx] >= ctg_supclusts->begs[supercluster_idx]) {
                    supercluster_idx++;
                }
                if (cluster2_idx >= int(calls2_vars->clusters.size()))
                    ERROR("Out of bounds cluster during write_results(): calls2")
                var2_idx++;
            }
        }
    }
    fclose(out_calls);
    
    //print truth variant information
    std::string out_truth_fn = g.out_prefix + "truth.tsv";
    INFO("  Printing truth variant results to '%s'", out_truth_fn.data());
    FILE* out_truth = fopen(out_truth_fn.data(), "w");
    fprintf(out_truth, "CONTIG\tPOS\tHAP\tREF\tALT\tQUAL\tTYPE\tERRTYPE\tCREDIT\tORIG_GT\tCLUSTER\tSUPERCLUSTER\tLOCATION\n");
    for (auto ctg_name : phasedata_ptr->contigs) {

        // set pointers to variants and superclusters
        std::shared_ptr<ctgVariants> truth1_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->truth1_vars;
        std::shared_ptr<ctgVariants> truth2_vars = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters->truth2_vars;
        std::shared_ptr<ctgClusters> ctg_supclusts = 
            phasedata_ptr->ctg_phasings[ctg_name].ctg_superclusters;

        int var1_idx = 0;
        int var2_idx = 0;
        int cluster1_idx = 0;
        int cluster2_idx = 0;
        int supercluster_idx = 0;
        while (var1_idx < truth1_vars->n || var2_idx < truth2_vars->n) {
            if (var2_idx >= truth2_vars->n || // only truth1 has remaining vars
                    truth1_vars->poss[var1_idx] < truth2_vars->poss[var2_idx]) { // truth1 var next
                if (truth1_vars->clusters[cluster1_idx] <= var1_idx) cluster1_idx++;
                fprintf(out_truth, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%s\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        truth1_vars->poss[var1_idx],
                        truth1_vars->haps[var1_idx],
                        truth1_vars->refs[var1_idx].data(),
                        truth1_vars->alts[var1_idx].data(),
                        truth1_vars->var_quals[var1_idx],
                        type_strs[truth1_vars->types[var1_idx]].data(),
                        error_strs[truth1_vars->errtypes[var1_idx]].data(),
                        truth1_vars->credit[var1_idx],
                        gt_strs[truth1_vars->orig_gts[var1_idx]].data(),
                        cluster1_idx,
                        supercluster_idx,
                        region_strs[truth1_vars->locs[var1_idx]].data()
                       );
                if (truth1_vars->poss[var1_idx] >= ctg_supclusts->begs[supercluster_idx]) {
                    supercluster_idx++;
                }
                if (cluster1_idx >= int(truth1_vars->clusters.size()))
                    ERROR("Out of bounds cluster during write_results(): truth1")
                var1_idx++;
            } else { // process truth2 var
                if (truth2_vars->clusters[cluster2_idx] <= var2_idx) cluster2_idx++;
                fprintf(out_truth, "%s\t%d\t%d\t%s\t%s\t%.2f\t%s\t%s\t%f\t%s\t%d\t%d\t%s\n",
                        ctg_name.data(),
                        truth2_vars->poss[var2_idx],
                        truth2_vars->haps[var2_idx],
                        truth2_vars->refs[var2_idx].data(),
                        truth2_vars->alts[var2_idx].data(),
                        truth2_vars->var_quals[var2_idx],
                        type_strs[truth2_vars->types[var2_idx]].data(),
                        error_strs[truth2_vars->errtypes[var2_idx]].data(),
                        truth2_vars->credit[var2_idx],
                        gt_strs[truth2_vars->orig_gts[var2_idx]].data(),
                        cluster2_idx,
                        supercluster_idx,
                        region_strs[truth2_vars->locs[var2_idx]].data()
                       );
                if (truth2_vars->poss[var2_idx] >= ctg_supclusts->begs[supercluster_idx]) {
                    supercluster_idx++;
                }
                if (cluster2_idx >= int(truth2_vars->clusters.size()))
                    ERROR("Out of bounds cluster during write_results(): truth2")
                var2_idx++;
            }
        }
    }
    fclose(out_truth);

    return;
}
