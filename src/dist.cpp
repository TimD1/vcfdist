#include <string>
#include <vector>
#include <cstdio>
#include <chrono>
#include <utility>
#include <queue>

#include "dist.h"
#include "print.h"


variantData edit_dist_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref) {

    // copy vcf header data over to results vcf
    variantData results;
    results.sample = vcf->sample;
    results.contigs = vcf->contigs;
    results.lengths = vcf->lengths;
    results.ref = vcf->ref;
    for (auto ctg : results.contigs) 
        for (int hap = 0; hap < 2; hap++) 
            results.ctg_variants[hap][ctg] = std::shared_ptr<ctgVariants>(new ctgVariants());

    // iterate over each haplotype
    int old_subs = 0;
    int old_inss = 0;
    int old_dels = 0;
    int old_edits = 0;
    int new_subs = 0;
    int new_inss = 0;
    int new_dels = 0;
    int new_edits = 0;
    for (int h = 0; h < 2; h++) {

        // iterate over each contig
        for (auto itr = vcf->ctg_variants[h].begin(); 
                itr != vcf->ctg_variants[h].end(); itr++) {
            std::string ctg = itr->first;
            std::shared_ptr<ctgVariants> vars = itr->second;
            if (vars->poss.size() == 0) continue;

            // iterate over each cluster of variants
            for (size_t cluster = 0; cluster < vars->clusters.size()-1; cluster++) {
                int beg_idx = vars->clusters[cluster];
                int end_idx = vars->clusters[cluster+1];
                int beg = vars->poss[beg_idx]-1;
                int end = vars->poss[end_idx-1] + vars->rlens[end_idx-1]+1;
                int old_subs_cluster = 0;
                int old_inss_cluster = 0;
                int old_dels_cluster = 0;

                // iterate over variants, summing edit distance
                for (int var = beg_idx; var < end_idx; var++) {
                    switch (vars->types[var]) {
                        case TYPE_SUB: old_subs_cluster++; break;
                        case TYPE_INS: old_inss_cluster += vars->alts[var].size(); break;
                        case TYPE_DEL: old_dels_cluster += vars->refs[var].size(); break;
                        default: ERROR("unexpected variant type (%i)", vars->types[var]) 
                                 std::exit(1); break;
                    }
                }

                // colored alignment
                int var = beg_idx;
                std::string ref_out_str = "";
                std::string alt_out_str = "";
                std::string alt_str = "";
                for (int ref_pos = beg; ref_pos < end;) {
                    if (ref_pos == vars->poss[var]) { // in variant
                        switch (vars->types[var]) {
                            case TYPE_INS:
                                alt_str += vars->alts[var];
                                alt_out_str += GREEN(vars->alts[var]);
                                ref_out_str += std::string(vars->alts[var].size(), ' ');
                                break;
                            case TYPE_DEL:
                                alt_out_str += std::string(vars->refs[var].size(), ' ');
                                ref_out_str += RED(vars->refs[var]);
                                ref_pos += vars->refs[var].size();
                                break;
                            case TYPE_SUB:
                                alt_str += vars->alts[var];
                                alt_out_str += " " + GREEN(vars->alts[var]);
                                ref_out_str += RED(vars->refs[var]) + " ";
                                ref_pos++;
                                break;
                            case TYPE_GRP:
                                alt_str += vars->alts[var];
                                alt_out_str += std::string(vars->refs[var].size(), ' ') 
                                    + GREEN(vars->alts[var]);
                                ref_out_str += RED(vars->refs[var]) + std::string(
                                        vars->alts[var].size(), ' ');
                                ref_pos += vars->refs[var].size();
                                break;
                        }
                        var++; // next variant
                    }
                    else { // match
                        try {
                            alt_str += ref->fasta.at(ctg)[ref_pos];
                            ref_out_str += ref->fasta.at(ctg)[ref_pos];
                            alt_out_str += ref->fasta.at(ctg)[ref_pos];
                            ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // do alignment
                std::string ref_str = ref->fasta.at(ctg).substr(beg, end-beg);
                int ref_len = ref_str.size();
                int alt_len = alt_str.size();
                int curr_score = 0;
                std::vector< std::vector<int> > ptrs(alt_len, std::vector<int>(ref_len));
                std::vector< std::vector<bool> > done(alt_len, std::vector<bool>(ref_len));
                std::queue< std::pair<int,int> > curr_q, next_q;

                // init
                ptrs[0][0] = PTR_DIAG;
                done[0][0] = true;

                // set first diag (wavefront)
                curr_q.push({0,0});
                for(int n = 1; n < std::min(ref_len, alt_len); n++) {
                    if (ref_str[n] == alt_str[n]) {
                        curr_q.push({n,n});
                        ptrs[n][n] = PTR_DIAG;
                        done[n][n] = true;
                    } else  {
                        break;
                    }
                }

                // continue alignment waves until fully aligned
                while (!done[alt_len-1][ref_len-1]) {

                    // expand to next wavefront
                    while (!curr_q.empty()) {

                        // get current cell
                        std::pair<int,int> cell = curr_q.front(); curr_q.pop();
                        int alt_idx = cell.first;
                        int ref_idx = cell.second;

                        // expand diagonal sub, then diagonally
                        if (alt_idx+1 < alt_len && ref_idx+1 < ref_len &&
                                !done[alt_idx+1][ref_idx+1]) {
                            next_q.push({alt_idx+1, ref_idx+1});
                            ptrs[alt_idx+1][ref_idx+1] |= PTR_SUB;
                            int off = 1;
                            while (alt_idx+1+off < alt_len && 
                                    ref_idx+1+off < ref_len && 
                                    !done[alt_idx+1+off][ref_idx+1+off] &&
                                    alt_str[alt_idx+1+off] == ref_str[ref_idx+1+off]) {
                                next_q.push({alt_idx+1+off, ref_idx+1+off});
                                ptrs[alt_idx+1+off][ref_idx+1+off] |= PTR_DIAG;
                                off++;
                            }
                        }

                        // expand down, then diagonally
                        if (alt_idx+1 < alt_len && !done[alt_idx+1][ref_idx]) {
                            next_q.push({alt_idx+1, ref_idx});
                            ptrs[alt_idx+1][ref_idx] |= PTR_UP;
                            int off = 1;
                            while (alt_idx+1+off < alt_len && 
                                    ref_idx+off < ref_len && 
                                    !done[alt_idx+1+off][ref_idx+off] &&
                                    alt_str[alt_idx+1+off] == ref_str[ref_idx+off]) {
                                next_q.push({alt_idx+1+off, ref_idx+off});
                                ptrs[alt_idx+1+off][ref_idx+off] |= PTR_DIAG;
                                off++;
                            }
                        }

                        // expand right, then diagonally
                        if (ref_idx+1 < ref_len && !done[alt_idx][ref_idx+1]) {
                            next_q.push({alt_idx, ref_idx+1});
                            ptrs[alt_idx][ref_idx+1] |= PTR_LEFT;
                            int off = 1;
                            while (alt_idx+off < alt_len && 
                                    ref_idx+1+off < ref_len && 
                                    !done[alt_idx+off][ref_idx+1+off] &&
                                    alt_str[alt_idx+off] == ref_str[ref_idx+1+off]) {
                                next_q.push({alt_idx+off, ref_idx+1+off});
                                ptrs[alt_idx+off][ref_idx+1+off] |= PTR_DIAG;
                                off++;
                            }
                        }
                    }

                    // current queue empty, transfer next over
                    while (!next_q.empty()) {

                        // get current cell
                        std::pair<int,int> cell = next_q.front(); next_q.pop();
                        int alt_idx = cell.first;
                        int ref_idx = cell.second;

                        if (!done[alt_idx][ref_idx]) {
                            // add to current wavefront, if new
                            curr_q.push(cell);

                            // mark cell as traversed
                            done[alt_idx][ref_idx] = true;
                        }
                    }
                    curr_score++;
                }

                // backtrack: count S/I/D, color path, set CIGAR
                int alt_idx = alt_len-1;
                int ref_idx = ref_len-1;
                int new_subs_cluster = 0;
                int new_inss_cluster = 0;
                int new_dels_cluster = 0;
                std::vector<int> cig(ref_len + alt_len);
                int cig_ptr = cig.size()-1;
                while (alt_idx >= 0 || ref_idx >= 0) {
                    ptrs[alt_idx][ref_idx] |= LEFT_PATH; // color print path
                    if (ptrs[alt_idx][ref_idx] & PTR_DIAG) {
                        ref_idx--; alt_idx--;
                        cig[cig_ptr--] = PTR_DIAG;
                        cig[cig_ptr--] = PTR_DIAG;
                    } else if (ptrs[alt_idx][ref_idx] & PTR_SUB) {
                        ref_idx--; alt_idx--; new_subs_cluster++;
                        cig[cig_ptr--] = PTR_SUB;
                        cig[cig_ptr--] = PTR_SUB;
                    } else if (ptrs[alt_idx][ref_idx] & PTR_LEFT) {
                        ref_idx--; new_dels_cluster++;
                        cig[cig_ptr--] = PTR_LEFT;
                    } else if (ptrs[alt_idx][ref_idx] & PTR_UP) {
                        alt_idx--; new_inss_cluster++;
                        cig[cig_ptr--] = PTR_UP;
                    } else {
                        ERROR("no pointer during backtrack at (%d,%d)", alt_idx, ref_idx);
                    }
                }

                // update totals
                int old_edits_cluster = old_subs_cluster + 
                    old_inss_cluster + old_dels_cluster;
                int new_edits_cluster = new_subs_cluster + 
                    new_inss_cluster + new_dels_cluster;
                old_subs += old_subs_cluster;
                new_subs += new_subs_cluster;
                old_inss += old_inss_cluster;
                new_inss += new_inss_cluster;
                old_dels += old_dels_cluster;
                new_dels += new_dels_cluster;
                old_edits += old_edits_cluster;
                new_edits += new_edits_cluster;

                // debug print alignment
                if (g.print_verbosity > 1 && old_edits_cluster != new_edits_cluster) {
                    printf("\n\nREF: %s\n", ref_out_str.data());
                    printf("ALT: %s\n", alt_out_str.data());
                    print_ptrs(ptrs, alt_str, ref_str);
                    printf("OLD: SUBs=%d  INSs=%d  DELs=%d  Total=%d\n", old_subs_cluster, 
                            old_inss_cluster, old_dels_cluster, old_edits_cluster);
                    printf("NEW: SUBs=%d  INSs=%d  DELs=%d  Total=%d\n", new_subs_cluster, 
                            new_inss_cluster, new_dels_cluster, new_edits_cluster);
                }


                // write new alignment to VCF
                ref_idx = 0, alt_idx = 0;
                for(size_t cig_idx = 0; cig_idx < cig.size();) {
                    int indel_len = 0;
                    switch (cig[cig_idx]) {

                        case PTR_DIAG: // no variant, update pointers
                            cig_idx += 2;
                            ref_idx++;
                            alt_idx++;
                            break;

                        case PTR_SUB: // substitution
                            cig_idx += 2;
                            results.ctg_variants[h][ctg]->add_var(beg+ref_idx, 1, h, 
                                    TYPE_SUB, std::string(1,ref_str[ref_idx]), 
                                    std::string(1,alt_str[alt_idx]), 60, 60);
                            ref_idx++;
                            alt_idx++;
                            break;

                        case PTR_LEFT: // deletion
                            cig_idx++; indel_len++;

                            // multi-base deletion
                            while (cig_idx < cig.size() && cig[cig_idx] == PTR_LEFT) {
                                cig_idx++; indel_len++;
                            }
                            results.ctg_variants[h][ctg]->add_var(beg+ref_idx,
                                    indel_len, h, TYPE_DEL, 
                                    ref_str.substr(ref_idx, indel_len),
                                    "", 60, 60);
                            ref_idx += indel_len;
                            break;

                        case PTR_UP: // insertion
                            cig_idx++; indel_len++;

                            // multi-base insertion
                            while (cig_idx < cig.size() && cig[cig_idx] == PTR_UP) {
                                cig_idx++; indel_len++;
                            }
                            results.ctg_variants[h][ctg]->add_var(beg+ref_idx,
                                    0, h, TYPE_INS, "", 
                                    alt_str.substr(alt_idx, indel_len), 60, 60);
                            alt_idx += indel_len;
                            break;
                    }
                }


            } // cluster
        } // contig
    } // hap
    INFO("Edit distance reduced from %d to %d.", old_edits, new_edits);
    INFO("OLD: SUBs=%d  INSs=%d  DELs=%d  Total=%d", old_subs, old_inss, old_dels, old_edits);
    INFO("NEW: SUBs=%d  INSs=%d  DELs=%d  Total=%d", new_subs, new_inss, new_dels, new_edits);

    return results;

}

/******************************************************************************/

std::shared_ptr<clusterData> edit_dist(
        std::unique_ptr<variantData> & calls, 
        std::unique_ptr<variantData> & truth, 
        std::shared_ptr<fastaData> ref) {

    // iterate over each contig
    int distance = 0;
    std::shared_ptr<clusterData> clusterdata_ptr(new clusterData(calls->contigs));
    for (std::string ctg : calls->contigs) {

        // initialize variant call info
        std::shared_ptr<ctgVariants> calls1_vars = calls->ctg_variants[0][ctg];
        std::shared_ptr<ctgVariants> calls2_vars = calls->ctg_variants[1][ctg];
        std::shared_ptr<ctgVariants> truth1_vars;
        std::shared_ptr<ctgVariants> truth2_vars;
        try {
            truth1_vars = truth->ctg_variants[0][ctg];
            truth2_vars = truth->ctg_variants[1][ctg];
        } catch (const std::exception & e) {
            ERROR("truth VCF does not contain contig '%s'", ctg.data());
        }
        clusterdata_ptr->ctg_superclusters[ctg]->set_variants(
                calls1_vars, calls2_vars, truth1_vars, truth2_vars);

        // for each cluster of variants (merge calls and truth haps)
        size_t calls1_clust_beg_idx = 0;
        size_t calls2_clust_beg_idx = 0;
        size_t truth1_clust_beg_idx = 0;
        size_t truth2_clust_beg_idx = 0;
        while (calls1_clust_beg_idx < calls1_vars->clusters.size()-1 || // last cluster added for end
               calls2_clust_beg_idx < calls2_vars->clusters.size()-1 ||
               truth1_clust_beg_idx < truth1_vars->clusters.size()-1 ||
               truth2_clust_beg_idx < truth2_vars->clusters.size()-1) {

            // init: empty supercluster
            size_t calls1_clust_end_idx = calls1_clust_beg_idx;
            size_t calls2_clust_end_idx = calls2_clust_beg_idx;
            size_t truth1_clust_end_idx = truth1_clust_beg_idx;
            size_t truth2_clust_end_idx = truth2_clust_beg_idx;

            int calls1_pos = (calls1_clust_beg_idx < calls1_vars->clusters.size()-1) ? 
                    calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int calls2_pos = (calls2_clust_beg_idx < calls2_vars->clusters.size()-1) ? 
                    calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int truth1_pos = (truth1_clust_beg_idx < truth1_vars->clusters.size()-1) ? 
                    truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int truth2_pos = (truth2_clust_beg_idx < truth2_vars->clusters.size()-1) ? 
                    truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();

            // initialize cluster merging with first to start
            int curr_end_pos = 0;
            if (calls1_pos <= truth1_pos && calls1_pos <= calls2_pos && calls1_pos <= truth2_pos) {
                calls1_clust_end_idx += 1;
                curr_end_pos = calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]-1] +
                        calls1_vars->rlens[calls1_vars->clusters[calls1_clust_end_idx]-1] + 1;
                calls1_pos = (calls1_clust_end_idx < calls1_vars->clusters.size()-1) ? 
                    calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (calls2_pos <= truth1_pos && calls2_pos <= calls1_pos && 
                    calls2_pos <= truth2_pos) {
                calls2_clust_end_idx += 1;
                curr_end_pos = calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]-1] +
                        calls2_vars->rlens[calls2_vars->clusters[calls2_clust_end_idx]-1] + 1;
                calls2_pos = (calls2_clust_end_idx < calls2_vars->clusters.size()-1) ? 
                    calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (truth1_pos <= calls1_pos && truth1_pos <= calls2_pos && 
                    truth1_pos <= truth2_pos) {
                truth1_clust_end_idx += 1;
                curr_end_pos = truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]-1] +
                        truth1_vars->rlens[truth1_vars->clusters[truth1_clust_end_idx]-1] + 1;
                truth1_pos = (truth1_clust_end_idx < truth1_vars->clusters.size()-1) ? 
                    truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else {
                truth2_clust_end_idx += 1;
                curr_end_pos = truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] +
                        truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1] + 1;
                truth2_pos = (truth2_clust_end_idx < truth2_vars->clusters.size()-1) ? 
                    truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            }


            // keep expanding cluster while possible
            bool just_merged = true;
            while (just_merged) {
                just_merged = false;
                while (truth1_pos < curr_end_pos + g.gap) {
                    truth1_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]-1] + 
                            truth1_vars->rlens[truth1_vars->clusters[truth1_clust_end_idx]-1] + 1);
                    truth1_pos = (truth1_clust_end_idx < truth1_vars->clusters.size()-1) ? 
                        truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (truth2_pos < curr_end_pos + g.gap) {
                    truth2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] + 
                            truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1] + 1);
                    truth2_pos = (truth2_clust_end_idx < truth2_vars->clusters.size()-1) ? 
                        truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (calls1_pos < curr_end_pos + g.gap) {
                    calls1_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos, 
                            calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]-1] + 
                            calls1_vars->rlens[calls1_vars->clusters[calls1_clust_end_idx]-1] + 1);
                    calls1_pos = (calls1_clust_end_idx < calls1_vars->clusters.size()-1) ? 
                        calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (calls2_pos < curr_end_pos + g.gap) {
                    calls2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]-1] + 
                            calls2_vars->rlens[calls2_vars->clusters[calls2_clust_end_idx]-1] + 1);
                    calls2_pos = (calls2_clust_end_idx < calls2_vars->clusters.size()-1) ? 
                        calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
            }

            // get supercluster start/end positions (allowing empty haps)
            int beg_pos = std::numeric_limits<int>::max();
            int end_pos = -1;
            if (calls1_clust_end_idx - calls1_clust_beg_idx) { // calls1 vars present
                beg_pos = std::min(beg_pos, 
                        calls1_vars->poss[calls1_vars->clusters[calls1_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        calls1_vars->poss[calls1_vars->clusters[calls1_clust_end_idx]-1] + 
                        calls1_vars->rlens[calls1_vars->clusters[calls1_clust_end_idx]-1]+1);
            }
            if (calls2_clust_end_idx - calls2_clust_beg_idx) { // calls2 vars present
                beg_pos = std::min(beg_pos, 
                        calls2_vars->poss[calls2_vars->clusters[calls2_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        calls2_vars->poss[calls2_vars->clusters[calls2_clust_end_idx]-1] + 
                        calls2_vars->rlens[calls2_vars->clusters[calls2_clust_end_idx]-1]+1);
            }
            if (truth1_clust_end_idx - truth1_clust_beg_idx) { // truth1 vars present
                beg_pos = std::min(beg_pos, 
                        truth1_vars->poss[truth1_vars->clusters[truth1_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        truth1_vars->poss[truth1_vars->clusters[truth1_clust_end_idx]-1] + 
                        truth1_vars->rlens[truth1_vars->clusters[truth1_clust_end_idx]-1]+1);
            }
            if (truth2_clust_end_idx - truth2_clust_beg_idx) { // truth2 vars present
                beg_pos = std::min(beg_pos, 
                        truth2_vars->poss[truth2_vars->clusters[truth2_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        truth2_vars->poss[truth2_vars->clusters[truth2_clust_end_idx]-1] + 
                        truth2_vars->rlens[truth2_vars->clusters[truth2_clust_end_idx]-1]+1);
            }

            // generate ref string
            std::string ref_str = ref->fasta.at(ctg).substr(beg_pos, end_pos-beg_pos);

            // generate calls1 and calls2 strings and pointers
            int calls1_var_idx = calls1_vars->clusters[calls1_clust_beg_idx];
            int calls2_var_idx = calls2_vars->clusters[calls2_clust_beg_idx];
            std::string calls1 = "", calls2 = "", calls1_str = "", calls2_str = ""; 
            std::vector<int> calls1_ptrs, calls2_ptrs;
            for (int calls1_ref_pos = beg_pos, calls2_ref_pos = beg_pos; calls1_ref_pos < end_pos || calls2_ref_pos < end_pos; ) {

                // CONSIDER TRUTH1 ONLY, PRIOR REFERENCE POSITION
                if (calls1_ref_pos < calls2_ref_pos) {
                    if (calls1_ref_pos == calls1_vars->poss[calls1_var_idx]) { // in calls1 variant
                        switch (calls1_vars->types[calls1_var_idx]) {
                            case TYPE_INS:
                                calls1 += calls1_vars->alts[calls1_var_idx];
                                calls1_ptrs.insert(calls1_ptrs.end(), calls1_vars->alts[calls1_var_idx].size(), -1);
                                calls1_str += GREEN(calls1_vars->alts[calls1_var_idx]);
                                break;
                            case TYPE_DEL:
                                calls1_str += std::string(calls1_vars->refs[calls1_var_idx].size(), ' ');
                                calls1_ref_pos += calls1_vars->refs[calls1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                calls1 += calls1_vars->alts[calls1_var_idx];
                                calls1_ptrs.insert(calls1_ptrs.end(), calls1_vars->alts[calls1_var_idx].size(), -1);
                                calls1_str += GREEN(calls1_vars->alts[calls1_var_idx]);
                                calls1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                calls1 += calls1_vars->alts[calls1_var_idx];
                                calls1_ptrs.insert(calls1_ptrs.end(), calls1_vars->alts[calls1_var_idx].size(), -1);
                                calls1_str += std::string(calls1_vars->refs[calls1_var_idx].size(), ' ') 
                                    + GREEN(calls1_vars->alts[calls1_var_idx]);
                                calls1_ref_pos += calls1_vars->refs[calls1_var_idx].size();
                                break;
                        }
                        calls1_var_idx++; // next variant
                    } else { // no calls1 variant, in calls2 variant
                        try {
                            calls1 += ref->fasta.at(ctg)[calls1_ref_pos];
                            calls1_ptrs.push_back(-1);
                            calls1_str += ref->fasta.at(ctg)[calls1_ref_pos];
                            calls1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // CONSIDER TRUTH2 ONLY, PRIOR REFERENCE POSITION
                else if (calls2_ref_pos < calls1_ref_pos) {
                    if (calls2_ref_pos == calls2_vars->poss[calls2_var_idx]) { // in calls2 variant
                        switch (calls2_vars->types[calls2_var_idx]) {
                            case TYPE_INS:
                                calls2 += calls2_vars->alts[calls2_var_idx];
                                calls2_ptrs.insert(calls2_ptrs.end(), calls2_vars->alts[calls2_var_idx].size(), -1);
                                calls2_str += GREEN(calls2_vars->alts[calls2_var_idx]);
                                break;
                            case TYPE_DEL:
                                calls2_str += std::string(calls2_vars->refs[calls2_var_idx].size(), ' ');
                                calls2_ref_pos += calls2_vars->refs[calls2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                calls2 += calls2_vars->alts[calls2_var_idx];
                                calls2_ptrs.insert(calls2_ptrs.end(), calls2_vars->alts[calls2_var_idx].size(), -1);
                                calls2_str += GREEN(calls2_vars->alts[calls2_var_idx]);
                                calls2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                calls2 += calls2_vars->alts[calls2_var_idx];
                                calls2_ptrs.insert(calls2_ptrs.end(), calls2_vars->alts[calls2_var_idx].size(), -1);
                                calls2_str += std::string(calls2_vars->refs[calls2_var_idx].size(), ' ') 
                                    + GREEN(calls2_vars->alts[calls2_var_idx]);
                                calls2_ref_pos += calls2_vars->refs[calls2_var_idx].size();
                                break;
                        }
                        calls2_var_idx++; // next variant
                    } else { // match
                        try {
                            calls2 += ref->fasta.at(ctg)[calls2_ref_pos];
                            calls2_ptrs.push_back(-1);
                            calls2_str += ref->fasta.at(ctg)[calls2_ref_pos];
                            calls2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // REFERENCE POSITIONS MATCH! POTENTIAL TRANSITIONS
                else {
                    bool calls1_var = false;
                    if (calls1_ref_pos == calls1_vars->poss[calls1_var_idx]) { // in calls1 variant
                        calls1_var = true;
                        switch (calls1_vars->types[calls1_var_idx]) {
                            case TYPE_INS:
                                calls1 += calls1_vars->alts[calls1_var_idx];
                                calls1_ptrs.insert(calls1_ptrs.end(), calls1_vars->alts[calls1_var_idx].size(), -1);
                                calls1_str += GREEN(calls1_vars->alts[calls1_var_idx]);
                                break;
                            case TYPE_DEL:
                                calls1_str += std::string(calls1_vars->refs[calls1_var_idx].size(), ' ');
                                calls1_ref_pos += calls1_vars->refs[calls1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                calls1 += calls1_vars->alts[calls1_var_idx];
                                calls1_ptrs.insert(calls1_ptrs.end(), calls1_vars->alts[calls1_var_idx].size(), -1);
                                calls1_str += GREEN(calls1_vars->alts[calls1_var_idx]);
                                calls1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                calls1 += calls1_vars->alts[calls1_var_idx];
                                calls1_ptrs.insert(calls1_ptrs.end(), calls1_vars->alts[calls1_var_idx].size(), -1);
                                calls1_str += std::string(calls1_vars->refs[calls1_var_idx].size(), ' ') 
                                    + GREEN(calls1_vars->alts[calls1_var_idx]);
                                calls1_ref_pos += calls1_vars->refs[calls1_var_idx].size();
                                break;
                        }
                        calls1_var_idx++; // next variant
                    } 

                    bool calls2_var = false;
                    if (calls2_ref_pos == calls2_vars->poss[calls2_var_idx]) { // in calls2 variant
                        calls2_var = true;
                        switch (calls2_vars->types[calls2_var_idx]) {
                            case TYPE_INS:
                                calls2 += calls2_vars->alts[calls2_var_idx];
                                calls2_ptrs.insert(calls2_ptrs.end(), calls2_vars->alts[calls2_var_idx].size(), -1);
                                calls2_str += GREEN(calls2_vars->alts[calls2_var_idx]);
                                break;
                            case TYPE_DEL:
                                calls2_str += std::string(calls2_vars->refs[calls2_var_idx].size(), ' ');
                                calls2_ref_pos += calls2_vars->refs[calls2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                calls2 += calls2_vars->alts[calls2_var_idx];
                                calls2_ptrs.insert(calls2_ptrs.end(), calls2_vars->alts[calls2_var_idx].size(), -1);
                                calls2_str += GREEN(calls2_vars->alts[calls2_var_idx]);
                                calls2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                calls2 += calls2_vars->alts[calls2_var_idx];
                                calls2_ptrs.insert(calls2_ptrs.end(), calls2_vars->alts[calls2_var_idx].size(), -1);
                                calls2_str += std::string(calls2_vars->refs[calls2_var_idx].size(), ' ') 
                                    + GREEN(calls2_vars->alts[calls2_var_idx]);
                                calls2_ref_pos += calls2_vars->refs[calls2_var_idx].size();
                                break;
                        }
                        calls2_var_idx++; // next variant

                    } 
                    
                    // ONE HAPLOTYPE WAS A VARIANT, INVALID POINTERS
                    if (!calls1_var && calls2_var) {
                        try {
                            calls1 += ref->fasta.at(ctg)[calls1_ref_pos];
                            calls1_ptrs.push_back(-1);
                            calls1_str += ref->fasta.at(ctg)[calls1_ref_pos];
                            calls1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                    if (calls1_var && !calls2_var) {
                        try {
                            calls2 += ref->fasta.at(ctg)[calls2_ref_pos];
                            calls2_ptrs.push_back(-1);
                            calls2_str += ref->fasta.at(ctg)[calls2_ref_pos];
                            calls2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }

                    // BOTH MATCH REFERENCE, ADD POINTERS
                    if (!calls1_var && !calls2_var) { // add pointers
                        try {
                            calls1_ptrs.push_back(calls2.size());
                            calls2_ptrs.push_back(calls1.size());
                            calls1 += ref->fasta.at(ctg)[calls1_ref_pos];
                            calls1_str += BLUE(ref->fasta.at(ctg)[calls1_ref_pos]);
                            calls1_ref_pos++;
                            calls2 += ref->fasta.at(ctg)[calls2_ref_pos];
                            calls2_str += BLUE(ref->fasta.at(ctg)[calls2_ref_pos]);
                            calls2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }
            }

            // generate truth1 and truth2 strings and pointers
            int truth1_var_idx = truth1_vars->clusters[truth1_clust_beg_idx];
            int truth2_var_idx = truth2_vars->clusters[truth2_clust_beg_idx];
            std::string truth1 = "", truth2 = "", truth1_str = "", truth2_str = ""; 
            std::vector<int> truth1_ptrs, truth2_ptrs;
            for (int truth1_ref_pos = beg_pos, truth2_ref_pos = beg_pos; truth1_ref_pos < end_pos || truth2_ref_pos < end_pos; ) {

                // CONSIDER TRUTH1 ONLY, PRIOR REFERENCE POSITION
                if (truth1_ref_pos < truth2_ref_pos) {
                    if (truth1_ref_pos == truth1_vars->poss[truth1_var_idx]) { // in truth1 variant
                        switch (truth1_vars->types[truth1_var_idx]) {
                            case TYPE_INS:
                                truth1 += truth1_vars->alts[truth1_var_idx];
                                truth1_ptrs.insert(truth1_ptrs.end(), truth1_vars->alts[truth1_var_idx].size(), -1);
                                truth1_str += GREEN(truth1_vars->alts[truth1_var_idx]);
                                break;
                            case TYPE_DEL:
                                truth1_str += std::string(truth1_vars->refs[truth1_var_idx].size(), ' ');
                                truth1_ref_pos += truth1_vars->refs[truth1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                truth1 += truth1_vars->alts[truth1_var_idx];
                                truth1_ptrs.insert(truth1_ptrs.end(), truth1_vars->alts[truth1_var_idx].size(), -1);
                                truth1_str += GREEN(truth1_vars->alts[truth1_var_idx]);
                                truth1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                truth1 += truth1_vars->alts[truth1_var_idx];
                                truth1_ptrs.insert(truth1_ptrs.end(), truth1_vars->alts[truth1_var_idx].size(), -1);
                                truth1_str += std::string(truth1_vars->refs[truth1_var_idx].size(), ' ') 
                                    + GREEN(truth1_vars->alts[truth1_var_idx]);
                                truth1_ref_pos += truth1_vars->refs[truth1_var_idx].size();
                                break;
                        }
                        truth1_var_idx++; // next variant
                    } else { // no truth1 variant, in truth2 variant
                        try {
                            truth1 += ref->fasta.at(ctg)[truth1_ref_pos];
                            truth1_ptrs.push_back(-1);
                            truth1_str += ref->fasta.at(ctg)[truth1_ref_pos];
                            truth1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // CONSIDER TRUTH2 ONLY, PRIOR REFERENCE POSITION
                else if (truth2_ref_pos < truth1_ref_pos) {
                    if (truth2_ref_pos == truth2_vars->poss[truth2_var_idx]) { // in truth2 variant
                        switch (truth2_vars->types[truth2_var_idx]) {
                            case TYPE_INS:
                                truth2 += truth2_vars->alts[truth2_var_idx];
                                truth2_ptrs.insert(truth2_ptrs.end(), truth2_vars->alts[truth2_var_idx].size(), -1);
                                truth2_str += GREEN(truth2_vars->alts[truth2_var_idx]);
                                break;
                            case TYPE_DEL:
                                truth2_str += std::string(truth2_vars->refs[truth2_var_idx].size(), ' ');
                                truth2_ref_pos += truth2_vars->refs[truth2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                truth2 += truth2_vars->alts[truth2_var_idx];
                                truth2_ptrs.insert(truth2_ptrs.end(), truth2_vars->alts[truth2_var_idx].size(), -1);
                                truth2_str += GREEN(truth2_vars->alts[truth2_var_idx]);
                                truth2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                truth2 += truth2_vars->alts[truth2_var_idx];
                                truth2_ptrs.insert(truth2_ptrs.end(), truth2_vars->alts[truth2_var_idx].size(), -1);
                                truth2_str += std::string(truth2_vars->refs[truth2_var_idx].size(), ' ') 
                                    + GREEN(truth2_vars->alts[truth2_var_idx]);
                                truth2_ref_pos += truth2_vars->refs[truth2_var_idx].size();
                                break;
                        }
                        truth2_var_idx++; // next variant
                    } else { // match
                        try {
                            truth2 += ref->fasta.at(ctg)[truth2_ref_pos];
                            truth2_ptrs.push_back(-1);
                            truth2_str += ref->fasta.at(ctg)[truth2_ref_pos];
                            truth2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // REFERENCE POSITIONS MATCH! POTENTIAL TRANSITIONS
                else {
                    bool truth1_var = false;
                    if (truth1_ref_pos == truth1_vars->poss[truth1_var_idx]) { // in truth1 variant
                        truth1_var = true;
                        switch (truth1_vars->types[truth1_var_idx]) {
                            case TYPE_INS:
                                truth1 += truth1_vars->alts[truth1_var_idx];
                                truth1_ptrs.insert(truth1_ptrs.end(), truth1_vars->alts[truth1_var_idx].size(), -1);
                                truth1_str += GREEN(truth1_vars->alts[truth1_var_idx]);
                                break;
                            case TYPE_DEL:
                                truth1_str += std::string(truth1_vars->refs[truth1_var_idx].size(), ' ');
                                truth1_ref_pos += truth1_vars->refs[truth1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                truth1 += truth1_vars->alts[truth1_var_idx];
                                truth1_ptrs.insert(truth1_ptrs.end(), truth1_vars->alts[truth1_var_idx].size(), -1);
                                truth1_str += GREEN(truth1_vars->alts[truth1_var_idx]);
                                truth1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                truth1 += truth1_vars->alts[truth1_var_idx];
                                truth1_ptrs.insert(truth1_ptrs.end(), truth1_vars->alts[truth1_var_idx].size(), -1);
                                truth1_str += std::string(truth1_vars->refs[truth1_var_idx].size(), ' ') 
                                    + GREEN(truth1_vars->alts[truth1_var_idx]);
                                truth1_ref_pos += truth1_vars->refs[truth1_var_idx].size();
                                break;
                        }
                        truth1_var_idx++; // next variant
                    } 

                    bool truth2_var = false;
                    if (truth2_ref_pos == truth2_vars->poss[truth2_var_idx]) { // in truth2 variant
                        truth2_var = true;
                        switch (truth2_vars->types[truth2_var_idx]) {
                            case TYPE_INS:
                                truth2 += truth2_vars->alts[truth2_var_idx];
                                truth2_ptrs.insert(truth2_ptrs.end(), truth2_vars->alts[truth2_var_idx].size(), -1);
                                truth2_str += GREEN(truth2_vars->alts[truth2_var_idx]);
                                break;
                            case TYPE_DEL:
                                truth2_str += std::string(truth2_vars->refs[truth2_var_idx].size(), ' ');
                                truth2_ref_pos += truth2_vars->refs[truth2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                truth2 += truth2_vars->alts[truth2_var_idx];
                                truth2_ptrs.insert(truth2_ptrs.end(), truth2_vars->alts[truth2_var_idx].size(), -1);
                                truth2_str += GREEN(truth2_vars->alts[truth2_var_idx]);
                                truth2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                truth2 += truth2_vars->alts[truth2_var_idx];
                                truth2_ptrs.insert(truth2_ptrs.end(), truth2_vars->alts[truth2_var_idx].size(), -1);
                                truth2_str += std::string(truth2_vars->refs[truth2_var_idx].size(), ' ') 
                                    + GREEN(truth2_vars->alts[truth2_var_idx]);
                                truth2_ref_pos += truth2_vars->refs[truth2_var_idx].size();
                                break;
                        }
                        truth2_var_idx++; // next variant

                    } 
                    
                    // ONE HAPLOTYPE WAS A VARIANT, INVALID POINTERS
                    if (!truth1_var && truth2_var) {
                        try {
                            truth1 += ref->fasta.at(ctg)[truth1_ref_pos];
                            truth1_ptrs.push_back(-1);
                            truth1_str += ref->fasta.at(ctg)[truth1_ref_pos];
                            truth1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                    if (truth1_var && !truth2_var) {
                        try {
                            truth2 += ref->fasta.at(ctg)[truth2_ref_pos];
                            truth2_ptrs.push_back(-1);
                            truth2_str += ref->fasta.at(ctg)[truth2_ref_pos];
                            truth2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }

                    // BOTH MATCH REFERENCE, ADD POINTERS
                    if (!truth1_var && !truth2_var) { // add pointers
                        try {
                            truth1_ptrs.push_back(truth2.size());
                            truth2_ptrs.push_back(truth1.size());
                            truth1 += ref->fasta.at(ctg)[truth1_ref_pos];
                            truth1_str += BLUE(ref->fasta.at(ctg)[truth1_ref_pos]);
                            truth1_ref_pos++;
                            truth2 += ref->fasta.at(ctg)[truth2_ref_pos];
                            truth2_str += BLUE(ref->fasta.at(ctg)[truth2_ref_pos]);
                            truth2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }
            }

            // ALIGNMENT
            std::vector<int> s(4);
            std::vector<std::string> calls {calls1, calls1, calls2, calls2};
            std::vector<std::string> truth {truth1, truth2, truth1, truth2};
            std::vector<int> calls_lens = 
                    {int(calls1.size()), int(calls1.size()), int(calls2.size()), int(calls2.size())};
            std::vector<int> truth_lens = 
                    {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};
            std::vector< std::vector< std::vector<int> > > offs(4), ptrs(4);

            // for each combination of calls and truth
            for(int i = 0; i < 4; i++) {

                int mat_len = calls_lens[i] + truth_lens[i] - 1;
                offs[i].push_back(std::vector<int>(mat_len,-2));
                offs[i][0][calls_lens[i]-1] = -1;
                ptrs[i].push_back(std::vector<int>(mat_len,PTR_NONE));
                bool done = false;
                while (true) {

                    // EXTEND WAVEFRONT
                    for (int d = 0; d < mat_len; d++) {
                        int off = offs[i][s[i]][d];
                        int diag = d + 1 - calls_lens[i];

                        // don't allow starting from untouched cells
                        if (off == -2) continue;

                        // check that it's within matrix
                        if (diag + off + 1 < 0) continue;
                        if (off > calls_lens[i] - 1) continue;
                        if (diag + off > truth_lens[i] - 1) continue;

                        // extend
                        while (off < calls_lens[i] - 1 && 
                               diag + off < truth_lens[i] - 1) {
                            if (calls[i][off+1] == truth[i][diag+off+1]) off++;
                            else break;
                        }
                        offs[i][s[i]][d] = off;

                        // finish if done
                        if (off == calls_lens[i] - 1 && 
                            off + diag == truth_lens[i] - 1)
                        { done = true; break; }

                    }
                    if (done) break;


                    // NEXT WAVEFRONT
                    // add wavefront, fill edge cells
                    offs[i].push_back(std::vector<int>(mat_len, -2));
                    // bottom left cells
                    if (s[i]+1 == calls_lens[i]-1)
                        offs[i][s[i]+1][0] = s[i]+1;
                    // top right cells
                    if (s[i]+1 == mat_len-1)
                        offs[i][s[i]+1][mat_len-1] = s[i]+1;

                    ptrs[i].push_back(std::vector<int>(mat_len));
                    ptrs[i][s[i]+1][0] = PTR_UP;
                    ptrs[i][s[i]+1][s[i]+1] = PTR_LEFT;

                    // central cells
                    for (int d = 1; d < mat_len-1; d++) {
                        int offleft = offs[i][s[i]][d-1];
                        int offtop  = (offs[i][s[i]][d+1] == -2) ? 
                            -2 : offs[i][s[i]][d+1]+1;
                        int offdiag = (offs[i][s[i]][d] == -2) ? 
                            -2 : offs[i][s[i]][d]+1;
                        if (offdiag >= offtop && offdiag >= offleft) {
                            offs[i][s[i]+1][d] = offdiag;
                            ptrs[i][s[i]+1][d] = PTR_SUB;
                        } else if (offleft >= offtop) {
                            offs[i][s[i]+1][d] = offleft;
                            ptrs[i][s[i]+1][d] = PTR_LEFT;
                        } else {
                            offs[i][s[i]+1][d] = offtop;
                            ptrs[i][s[i]+1][d] = PTR_UP;
                        }
                    }
                    ++s[i];
                }
            }


            // PRINT RESULTS
            int dist = std::min(s[CALLS1_TRUTH1]+s[CALLS2_TRUTH2], 
                    s[CALLS2_TRUTH1]+s[CALLS1_TRUTH2]);
            if (g.print_verbosity >= 1 && dist) {
                // print cluster info
                printf("\n\nCALLS1: %zu clusters\n", calls1_clust_end_idx-calls1_clust_beg_idx);
                for(size_t i = calls1_clust_beg_idx; i < calls1_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, calls1_vars->clusters[i+1]-calls1_vars->clusters[i]);
                    for(int j = calls1_vars->clusters[i]; j < calls1_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), calls1_vars->poss[j], 
                                calls1_vars->refs[j].size() ? calls1_vars->refs[j].data() : "_", 
                                calls1_vars->alts[j].size() ? calls1_vars->alts[j].data() : "_");
                    }
                }
                printf("CALLS2: %zu clusters\n", calls2_clust_end_idx-calls2_clust_beg_idx);
                for(size_t i = calls2_clust_beg_idx; i < calls2_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, calls2_vars->clusters[i+1]-calls2_vars->clusters[i]);
                    for(int j = calls2_vars->clusters[i]; j < calls2_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), calls2_vars->poss[j], 
                                calls2_vars->refs[j].size() ? calls2_vars->refs[j].data() : "_", 
                                calls2_vars->alts[j].size() ? calls2_vars->alts[j].data() : "_");
                    }
                }
                printf("TRUTH1: %zu clusters\n", truth1_clust_end_idx-truth1_clust_beg_idx);
                for(size_t i = truth1_clust_beg_idx; i < truth1_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, truth1_vars->clusters[i+1]-truth1_vars->clusters[i]);
                    for(int j = truth1_vars->clusters[i]; j < truth1_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), truth1_vars->poss[j], 
                                truth1_vars->refs[j].size() ? truth1_vars->refs[j].data() : "_", 
                                truth1_vars->alts[j].size() ? truth1_vars->alts[j].data() : "_");
                    }
                }
                printf("TRUTH2: %zu clusters\n", truth2_clust_end_idx-truth2_clust_beg_idx);
                for(size_t i = truth2_clust_beg_idx; i < truth2_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, truth2_vars->clusters[i+1]-truth2_vars->clusters[i]);
                    for(int j = truth2_vars->clusters[i]; j < truth2_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), truth2_vars->poss[j], 
                                truth2_vars->refs[j].size() ? truth2_vars->refs[j].data() : "_", 
                                truth2_vars->alts[j].size() ? truth2_vars->alts[j].data() : "_");
                    }
                }
                
                printf("ORIG: %s\n", ref_str.data());
                printf("CALLS1: %s\n", calls1_str.data());
                printf("CALLS2: %s\n", calls2_str.data());
                printf("TRUTH1: %s\n", truth1_str.data());
                printf("TRUTH2: %s\n", truth2_str.data());
                printf("Edit Distance: %d\n", dist);
            }

            // DEBUG PRINT
            if (g.print_verbosity >= 2 && dist) {
                for(int h = 0; h < 4; h++) { // 4 alignments
                    printf("\n%s ALIGNMENT (distance %d)\n", 
                            aln_strs[h].data(), s[h]);

                    // create array
                    std::vector< std::vector<char> > ptr_str;
                    for (int i = 0; i < calls_lens[h]; i++)
                        ptr_str.push_back(std::vector<char>(truth_lens[h], '.'));

                    // modify array with pointers
                    int mat_len = calls_lens[h] + truth_lens[h] - 1;
                    for (int si = 0; si <= s[h]; si++) {
                        for(int di = 0; di < mat_len; di++) {
                            int diag = di + 1 - calls_lens[h];
                            int off = offs[h][si][di];

                            // check that indices are within bounds
                            int calls_pos = off;
                            int truth_pos = diag + off;
                            if (calls_pos < 0 || truth_pos < 0) continue;
                            if (calls_pos > calls_lens[h]-1 || 
                                    truth_pos > truth_lens[h]-1) continue;

                            // special case: main diag, no previous edits
                            if (si == 0 && diag == 0) {
                                while (calls_pos >= 0) {
                                    ptr_str[calls_pos--][truth_pos--] = '\\';
                                }
                            }
                            // left edge
                            else if (calls_pos == 0) {
                                ptr_str[calls_pos][truth_pos] = '-';
                            } 
                            // top edge
                            else if (truth_pos == 0) {
                                ptr_str[calls_pos][truth_pos] = '|';
                            } 
                            else {
                                // follow diagonal
                                int top_off = (di < mat_len-1) ? offs[h][si-1][di+1]+1 : -2;
                                int left_off = (di > 0) ? offs[h][si-1][di-1] : -2;
                                int diag_off = offs[h][si-1][di]+1;
                                while (calls_pos > 0 && truth_pos > 0 && 
                                        calls_pos > top_off && 
                                        calls_pos > left_off && 
                                        calls_pos > diag_off) {
                                    ptr_str[calls_pos--][truth_pos--] = '\\';
                                }
                                // check left/up
                                if (calls_pos == diag_off) {
                                    ptr_str[calls_pos][truth_pos] = 'X';
                                }
                                else if (calls_pos == top_off) {
                                    ptr_str[calls_pos][truth_pos] = '|';
                                } else if (calls_pos == left_off) {
                                    ptr_str[calls_pos][truth_pos] = '-';
                                }
                            }
                        }
                    }

                    // print array
                    for (int i = -1; i < calls_lens[h]; i++) {
                        for (int j = -1; j < truth_lens[h]; j++) {
                            if (i < 0 && j < 0) {
                                printf("  ");
                            }
                            else if (i < 0) {
                                printf("%c", truth[h][j]);
                            } else if (j < 0) {
                                printf("\n%c ", calls[h][i]);
                            } else {
                                printf("%c", ptr_str[i][j]);
                            }
                        }
                    }
                    printf("\n");

                } // 4 alignments
            } // debug print

            // get cluster phasing
            int orig_phase_dist = s[CALLS1_TRUTH1] + s[CALLS2_TRUTH2];
            int swap_phase_dist = s[CALLS2_TRUTH1] + s[CALLS1_TRUTH2];
            int phase = PHASE_NONE; // default either way if equal dist
            if (orig_phase_dist < swap_phase_dist) phase = PHASE_ORIG;
            if (swap_phase_dist < orig_phase_dist) phase = PHASE_SWAP;

            // save alignment information
            clusterdata_ptr->ctg_superclusters[ctg]->add_supercluster(
                    calls1_clust_beg_idx, calls1_clust_end_idx,
                    calls2_clust_beg_idx, calls2_clust_end_idx,
                    truth1_clust_beg_idx, truth1_clust_end_idx,
                    truth2_clust_beg_idx, truth2_clust_end_idx,
                    beg_pos, end_pos,
                    phase, orig_phase_dist, swap_phase_dist);

            // reset for next merged cluster
            calls1_clust_beg_idx = calls1_clust_end_idx;
            calls2_clust_beg_idx = calls2_clust_end_idx;
            truth1_clust_beg_idx = truth1_clust_end_idx;
            truth2_clust_beg_idx = truth2_clust_end_idx;
            distance += std::min(orig_phase_dist, swap_phase_dist);

        } // each cluster
    } // each contig
    INFO("Total edit distance: %d", distance);

    return clusterdata_ptr;
}
