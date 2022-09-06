#include <string>
#include <vector>
#include <cstdio>
#include <chrono>
#include <utility>
#include <queue>

#include "dist.h"
#include "print.h"
#include "cluster.h"


vcfData edit_dist_realign(const vcfData* vcf, const fastaData* const ref) {

    // copy vcf header data over to results vcf
    vcfData results;
    results.sample = vcf->sample;
    results.contigs = vcf->contigs;
    results.lengths = vcf->lengths;
    results.ref = vcf->ref;
    for (auto ctg : results.contigs) 
        for (int hap = 0; hap < 2; hap++) 
            results.hapcalls[hap][ctg] = variantCalls();

    // iterate over each haplotype
    int clusters = 0;
    int new_ed_clusters = 0;
    int old_ed = 0;
    int new_ed = 0;
    for (int h = 0; h < 2; h++) {

        // iterate over each contig
        for (auto itr = vcf->hapcalls[h].begin(); 
                itr != vcf->hapcalls[h].end(); itr++) {
            std::string ctg = itr->first;
            variantCalls vars = itr->second;
            if (vars.poss.size() == 0) continue;

            // iterate over each cluster of variants
            for (size_t cluster = 0; cluster < vars.clusters.size()-1; cluster++) {
                int beg_idx = vars.clusters[cluster];
                int end_idx = vars.clusters[cluster+1];
                int beg = vars.poss[beg_idx]-1;
                int end = vars.poss[end_idx-1] + vars.rlens[end_idx-1]+1;
                int subs = 0;
                int inss = 0;
                int dels = 0;

                // iterate over variants, summing edit distance
                for (int var = beg_idx; var < end_idx; var++) {
                    switch (vars.types[var]) {
                        case TYPE_SUB: subs++; break;
                        case TYPE_INS: inss += vars.alts[var].size(); break;
                        case TYPE_DEL: dels += vars.refs[var].size(); break;
                        default: ERROR("unexpected variant type (%i)", vars.types[var]) 
                                 std::exit(1); break;
                    }
                }

                // colored alignment
                int var = beg_idx;
                std::string ref_out_str = "";
                std::string alt_out_str = "";
                std::string alt_str = "";
                for (int ref_pos = beg; ref_pos < end;) {
                    if (ref_pos == vars.poss[var]) { // in variant
                        switch (vars.types[var]) {
                            case TYPE_INS:
                                alt_str += vars.alts[var];
                                alt_out_str += GREEN(vars.alts[var]);
                                ref_out_str += std::string(vars.alts[var].size(), ' ');
                                break;
                            case TYPE_DEL:
                                alt_out_str += std::string(vars.refs[var].size(), ' ');
                                ref_out_str += RED(vars.refs[var]);
                                ref_pos += vars.refs[var].size();
                                break;
                            case TYPE_SUB:
                                alt_str += vars.alts[var];
                                alt_out_str += " " + GREEN(vars.alts[var]);
                                ref_out_str += RED(vars.refs[var]) + " ";
                                ref_pos++;
                                break;
                            case TYPE_GRP:
                                alt_str += vars.alts[var];
                                alt_out_str += std::string(vars.refs[var].size(), ' ') 
                                    + GREEN(vars.alts[var]);
                                ref_out_str += RED(vars.refs[var]) + std::string(
                                        vars.alts[var].size(), ' ');
                                ref_pos += vars.refs[var].size();
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

                // FORWARD PASS ALIGNMENT
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

                // left backtrack (prioritize LEFT)
                int alt_idx = alt_len-1;
                int ref_idx = ref_len-1;
                std::vector< std::pair<int,int> > left_path;
                while (alt_idx >= 0 || ref_idx >= 0) {
                    left_path.push_back({alt_idx, ref_idx});
                    ptrs[alt_idx][ref_idx] |= LEFT_PATH;
                    if (ptrs[alt_idx][ref_idx] & PTR_LEFT) {
                        ref_idx--;
                    } else if (ptrs[alt_idx][ref_idx] & PTR_DIAG) {
                        ref_idx--; alt_idx--;
                    } else if (ptrs[alt_idx][ref_idx] & PTR_UP) {
                        alt_idx--;
                    } else {
                        ERROR("no pointer during left backtrack at (%d,%d)", alt_idx, ref_idx);
                    }
                }
                std::reverse(left_path.begin(), left_path.end());

                // right backtrack (prioritize UP)
                alt_idx = alt_len-1;
                ref_idx = ref_len-1;
                std::vector< std::pair<int,int> > right_path;
                while (alt_idx >= 0 || ref_idx >= 0) {
                    right_path.push_back({alt_idx, ref_idx});
                    ptrs[alt_idx][ref_idx] |= RIGHT_PATH;
                    if (ptrs[alt_idx][ref_idx] & PTR_UP) {
                        alt_idx--;
                    } else if (ptrs[alt_idx][ref_idx] & PTR_DIAG) {
                        ref_idx--; alt_idx--;
                    } else if (ptrs[alt_idx][ref_idx] & PTR_LEFT) {
                        ref_idx--;
                    } else {
                        ERROR("no pointer during right backtrack at (%d,%d)", alt_idx, ref_idx);
                    }
                }
                std::reverse(right_path.begin(), right_path.end());
                

                // BACKWARD PASS ALIGNMENT (MIN EDIT EVENTS)
                int curr_edits = 0;
                std::vector< std::vector<int> > bt_ptrs(alt_len, std::vector<int>(ref_len));
                curr_q = {}; next_q = {};
                std::queue< std::pair<int,int> > sub_q;

                // set first diag (wavefront)
                alt_idx = alt_len-1;
                ref_idx = ref_len-1;
                while (alt_idx >= 0 && ref_idx >= 0 && 
                        ptrs[alt_idx][ref_idx] & PTR_DIAG) {
                    bt_ptrs[alt_idx][ref_idx] |= (PTR_DIAG | PTR_DONE);
                    alt_idx--; ref_idx--;
                }
                if (alt_idx >= 0 && ref_idx >= 0) {
                    bt_ptrs[alt_idx][ref_idx] |= PTR_DONE;
                    curr_q.push({alt_idx, ref_idx});
                    sub_q.push({alt_idx, ref_idx});
                }

                // continue alignment waves until fully aligned
                while (!(bt_ptrs[0][0] & PTR_DONE)) {

                    // expand to next wavefront
                    while (!curr_q.empty()) {

                        // get current cell
                        std::pair<int,int> cell = curr_q.front(); curr_q.pop();
                        alt_idx = cell.first;
                        ref_idx = cell.second;

                        // expand up, then diagonally
                        int alt_idx_ins = alt_idx;
                        int ref_idx_ins = ref_idx;
                        while (alt_idx_ins > 0 && 
                                !(bt_ptrs[alt_idx_ins-1][ref_idx_ins] & PTR_DONE) &&
                                ptrs[alt_idx_ins][ref_idx_ins] & PTR_UP) {
                            bt_ptrs[alt_idx_ins][ref_idx_ins] |= PTR_UP;
                            int alt_idx_ins_slide = alt_idx_ins - 1;
                            int ref_idx_ins_slide = ref_idx_ins;
                            while (alt_idx_ins_slide > 0 && ref_idx_ins_slide > 0 &&
                                    !(bt_ptrs[alt_idx_ins_slide-1][ref_idx_ins_slide-1] & PTR_DONE) &&
                                    !(bt_ptrs[alt_idx_ins_slide-1][ref_idx_ins_slide-1] & PTR_NEXT && 
                                        ptrs[alt_idx_ins_slide-1][ref_idx_ins_slide-1] & PTR_UP && 
                                        ptrs[alt_idx_ins_slide][ref_idx_ins_slide-1] & PTR_UP) &&
                                    !(bt_ptrs[alt_idx_ins_slide-1][ref_idx_ins_slide-1] & PTR_NEXT && 
                                        ptrs[alt_idx_ins_slide-1][ref_idx_ins_slide-1] & PTR_LEFT && 
                                        ptrs[alt_idx_ins_slide-1][ref_idx_ins_slide] & PTR_LEFT) &&
                                    ptrs[alt_idx_ins_slide][ref_idx_ins_slide] & PTR_DIAG) {
                                bt_ptrs[alt_idx_ins_slide][ref_idx_ins_slide] |= (PTR_DIAG | PTR_NEXT);
                                alt_idx_ins_slide--; ref_idx_ins_slide--;
                            }
                            if (alt_idx_ins_slide >= 0 && ref_idx_ins_slide >= 0) {
                                bt_ptrs[alt_idx_ins_slide][ref_idx_ins_slide] |= PTR_NEXT;
                                next_q.push({alt_idx_ins_slide, ref_idx_ins_slide});
                            }
                            alt_idx_ins--;
                        }

                        // expand left, then diagonally
                        int alt_idx_del = alt_idx;
                        int ref_idx_del = ref_idx;
                        while (ref_idx_del > 0 && 
                                !(bt_ptrs[alt_idx_del][ref_idx_del-1] & PTR_DONE) &&
                                ptrs[alt_idx_del][ref_idx_del] & PTR_LEFT) {
                            bt_ptrs[alt_idx_del][ref_idx_del] |= PTR_LEFT;
                            int alt_idx_del_slide = alt_idx_del;
                            int ref_idx_del_slide = ref_idx_del - 1;
                            while (alt_idx_del_slide > 0 && ref_idx_del_slide > 0 &&
                                    !(bt_ptrs[alt_idx_del_slide-1][ref_idx_del_slide-1] & PTR_DONE) &&
                                    !(bt_ptrs[alt_idx_del_slide-1][ref_idx_del_slide-1] & PTR_NEXT && 
                                        ptrs[alt_idx_del_slide-1][ref_idx_del_slide-1] & PTR_UP && 
                                        ptrs[alt_idx_del_slide][ref_idx_del_slide-1] & PTR_UP) &&
                                    !(bt_ptrs[alt_idx_del_slide-1][ref_idx_del_slide-1] & PTR_NEXT && 
                                        ptrs[alt_idx_del_slide-1][ref_idx_del_slide-1] & PTR_LEFT && 
                                        ptrs[alt_idx_del_slide-1][ref_idx_del_slide] & PTR_LEFT) &&
                                        ptrs[alt_idx_del_slide][ref_idx_del_slide] & PTR_DIAG) {
                                bt_ptrs[alt_idx_del_slide][ref_idx_del_slide] |= (PTR_DIAG | PTR_NEXT);
                                alt_idx_del_slide--; ref_idx_del_slide--;
                            }
                            if (alt_idx_del_slide >= 0 && ref_idx_del_slide >= 0) {
                                bt_ptrs[alt_idx_del_slide][ref_idx_del_slide] |= PTR_NEXT;
                                next_q.push({alt_idx_del_slide, ref_idx_del_slide});
                            }
                            ref_idx_del--;
                        }
                    }

                    while (!sub_q.empty()) {
                        // get current cell
                        std::pair<int,int> cell = sub_q.front(); sub_q.pop();
                        alt_idx = cell.first;
                        ref_idx = cell.second;

                        // expand sub, then continue diagonally
                        if (alt_idx > 0 && ref_idx > 0 && 
                                !(ptrs[alt_idx][ref_idx] & PTR_DIAG) &&
                                !(bt_ptrs[alt_idx-1][ref_idx-1] & PTR_DONE) &&
                                !(bt_ptrs[alt_idx-1][ref_idx-1] & PTR_NEXT && 
                                    ptrs[alt_idx-1][ref_idx-1] & PTR_UP && 
                                    ptrs[alt_idx][ref_idx-1] & PTR_UP) &&
                                !(bt_ptrs[alt_idx-1][ref_idx-1] & PTR_NEXT && 
                                    ptrs[alt_idx-1][ref_idx-1] & PTR_LEFT && 
                                    ptrs[alt_idx-1][ref_idx] & PTR_LEFT)
                                ) {
                            bt_ptrs[alt_idx][ref_idx] |= PTR_SUB;
                            int alt_idx_sub = alt_idx - 1;
                            int ref_idx_sub = ref_idx - 1;
                            while (alt_idx_sub > 0 && ref_idx_sub > 0 &&
                                    !(bt_ptrs[alt_idx_sub-1][ref_idx_sub-1] & PTR_DONE) &&
                                    ptrs[alt_idx_sub][ref_idx_sub] & PTR_DIAG) {
                                bt_ptrs[alt_idx_sub][ref_idx_sub] |= (PTR_DIAG | PTR_NEXT);
                                alt_idx_sub--; ref_idx_sub--;
                            }
                            if (alt_idx_sub >= 0 && ref_idx_sub >= 0) {
                                bt_ptrs[alt_idx_sub][ref_idx_sub] |= PTR_NEXT;
                                next_q.push({alt_idx_sub, ref_idx_sub});
                            }
                        }

                    }

                    /* if (cluster == 43531) { */
                    /*     print_ptrs(ptrs, alt_str, ref_str); */
                    /*     print_ptrs(bt_ptrs, alt_str, ref_str); */
                    /* } */

                    // current queue empty, transfer next over
                    while (!next_q.empty()) {

                        // get current cell
                        std::pair<int,int> cell = next_q.front(); next_q.pop();
                        int alt_idx = cell.first;
                        int ref_idx = cell.second;

                        if (alt_idx >= 0 && ref_idx >= 0 && 
                                !(bt_ptrs[alt_idx][ref_idx] & PTR_DONE)) {
                            bt_ptrs[alt_idx][ref_idx] |= PTR_DONE;
                            curr_q.push(cell);
                            sub_q.push(cell);
                        } 
                    }

                    // update completed cells
                    for (alt_idx = 0; alt_idx < alt_len; alt_idx++) {
                        for (ref_idx = 0; ref_idx < ref_len; ref_idx++) {
                            if (bt_ptrs[alt_idx][ref_idx] & PTR_NEXT) {
                                bt_ptrs[alt_idx][ref_idx] |= PTR_DONE;
                                bt_ptrs[alt_idx][ref_idx] &= (~PTR_NEXT);
                            }
                        }
                    }

                    curr_edits++;
                }
                bt_ptrs[0][0] |= PTR_DIAG;

                /* if (g.print_verbosity >= 2 && 2*subs + inss + dels != curr_score) { */
                /*     print_ptrs(ptrs, alt_str, ref_str); */
                /*     print_ptrs(bt_ptrs, alt_str, ref_str); */
                /* } */

                // forward-pass finish
                ref_idx = 0; alt_idx = 0;
                int new_subs = 0;
                int new_inss = 0;
                int new_dels = 0;
                while (ref_idx < ref_len || alt_idx < alt_len) {
                    if (bt_ptrs[alt_idx][ref_idx] & PTR_DIAG) {
                        bt_ptrs[alt_idx][ref_idx] |= LEFT_PATH;
                        alt_idx++; ref_idx++;
                    } else if (bt_ptrs[alt_idx][ref_idx] & PTR_SUB) {
                        bt_ptrs[alt_idx][ref_idx] |= LEFT_PATH;
                        new_subs++; alt_idx++; ref_idx++;
                    } else if (ref_idx > 0 && bt_ptrs[alt_idx][ref_idx-1] & PTR_UP) {
                        bt_ptrs[alt_idx][ref_idx-1] |= LEFT_PATH;
                        new_inss++; alt_idx++;
                    } else if (alt_idx > 0 && bt_ptrs[alt_idx-1][ref_idx] & PTR_LEFT) {
                        bt_ptrs[alt_idx-1][ref_idx] |= LEFT_PATH;
                        new_dels++; ref_idx++;
                    } else {
                        ERROR("error during second forward pass");
                    }
                }

                // DEBUG PRINT EDIT FINDING
                if (g.print_verbosity >= 2 && subs + inss + dels != new_subs + new_inss + new_dels) 
                /* if (g.print_verbosity >= 2 && 2*subs + inss + dels != curr_score) */ 
                {
                    std::string var_ref;
                    std::string var_alt;
                    printf("\n\n  Cluster %i: %s:%d-%d\n\n"
                            "    old edit distance: %d, %d variants (%dS %dI %dD)\n",
                            int(cluster), ctg.data(), beg, end, 
                            subs*2 + inss + dels, end_idx-beg_idx, subs, inss, dels);
                    for (int var = beg_idx; var < end_idx; var++) {
                        if (vars.refs[var].size() == 0) var_ref = "-"; 
                        else var_ref = vars.refs[var];
                        if (vars.alts[var].size() == 0) var_alt = "-"; 
                        else var_alt = vars.alts[var];
                        printf("      %s:%i hap%i %s\t%s\t%s\n", 
                                ctg.data(), vars.poss[var], h, 
                                type_strs[vars.types[var]].data(), 
                                var_ref.data(), var_alt.data());
                    }
                    printf("      REF: %s\n      ALT: %s\n", 
                            ref_out_str.data(), alt_out_str.data());
                    printf("\n    new edit distance: %d, %d variants (%dS %dI %dD)\n", 
                            curr_score, curr_edits, new_subs, new_inss, new_dels);
                    print_ptrs(ptrs, alt_str, ref_str);
                    print_ptrs(bt_ptrs, alt_str, ref_str);
                }

                // init (first base always matches)
                size_t left_idx = 1;
                size_t leftmost_idx = 1;
                size_t right_idx = 1;

                // follow both entire paths
                while (left_idx < left_path.size() || right_idx < right_path.size()) {

                    // follow matching paths (we know first base matches)
                    while (left_path[left_idx] == right_path[right_idx] && 
                            left_path[left_idx].first - left_path[left_idx-1].first == 1 &&
                            left_path[left_idx].second - left_path[left_idx-1].second == 1) {
                        left_idx++; right_idx++;
                    }
                    if (left_idx >= left_path.size() && right_idx >= right_path.size()) break;

                    // init stats for diverged region
                    int area = 0;
                    int grp_inss = 0;
                    int grp_dels = 0;

                    // get right-most position of left path in this row
                    leftmost_idx = left_idx;
                    while (left_idx < left_path.size() - 1 && 
                            left_path[left_idx+1].first == left_path[left_idx].first) {

                        // count inss/dels on this path
                        if (left_path[left_idx].first > left_path[left_idx-1].first &&
                            left_path[left_idx].second == left_path[left_idx-1].second)
                            grp_inss++;
                        if (left_path[left_idx].second > left_path[left_idx-1].second &&
                            left_path[left_idx].first == left_path[left_idx-1].first)
                            grp_dels++;

                        left_idx++;
                    }
                    if (left_path[left_idx].first > left_path[left_idx-1].first &&
                        left_path[left_idx].second == left_path[left_idx-1].second)
                        grp_inss++;
                    if (left_path[left_idx].second > left_path[left_idx-1].second &&
                        left_path[left_idx].first == left_path[left_idx-1].first)
                        grp_dels++;

                    // get start of diverged region
                    int alt_beg = left_path[leftmost_idx-1].first;
                    int ref_beg = left_path[leftmost_idx-1].second;

                    // path has diverged; continue until merged
                    bool merged = false;
                    while (!merged) {

                        // follow right path until same row as left
                        while(left_path[left_idx].first > right_path[right_idx].first)
                            right_idx++;
                        area += std::max(0, right_path[right_idx].second - 
                                left_path[leftmost_idx].second);

                        // check if paths merge
                        if (left_path[left_idx] == right_path[right_idx]) { // next cell same
                            if ((left_idx == left_path.size() - 1 &&  // matrix end
                                right_idx == right_path.size() - 1) || // next step diagonal for l/r
                                (left_path[left_idx+1] == right_path[right_idx+1] &&
                                 left_path[left_idx+1].first - left_path[left_idx].first == 1 &&
                                 left_path[left_idx+1].second - left_path[left_idx].second == 1)) {
                                merged = true;

                                /* int alt_end = left_path[left_idx].first + 1; */
                                /* int ref_end = left_path[left_idx].second + 1; */
                                /* if (subs*2 + inss + dels != curr_score) */ 
                                /* printf("I=%d D=%d ref[%d:%d]=%s alt[%d:%d]=%s, area=%d\n", */
                                /*         grp_inss, grp_dels, ref_beg, ref_end, */ 
                                /*         ref_str.substr(ref_beg, ref_end-ref_beg).data(), */
                                /*         alt_beg, alt_end, */ 
                                /*         alt_str.substr(alt_beg, alt_end-alt_beg).data(), area); */

                            }
                        }

                        // no merge, next positions
                        right_idx++;
                        left_idx++; 

                        // get right-most position of left path in this row
                        leftmost_idx = left_idx;
                        while (left_idx < left_path.size() - 1 && 
                                left_path[left_idx+1].first == left_path[left_idx].first) {

                            // count inss/dels on this path
                            if (left_path[left_idx].first > left_path[left_idx-1].first &&
                                left_path[left_idx].second == left_path[left_idx-1].second)
                                grp_inss++;
                            if (left_path[left_idx].second > left_path[left_idx-1].second &&
                                left_path[left_idx].first == left_path[left_idx-1].first)
                                grp_dels++;

                            left_idx++;
                        }
                        if (left_path[left_idx].first > left_path[left_idx-1].first &&
                            left_path[left_idx].second == left_path[left_idx-1].second)
                            grp_inss++;
                        if (left_path[left_idx].second > left_path[left_idx-1].second &&
                            left_path[left_idx].first == left_path[left_idx-1].first)
                            grp_dels++;
                    }
                }

                /* // update counters */
                /* old_ed += subs*2 + inss + dels; */
                /* new_ed += curr_score; */
                /* if (subs*2 + inss + dels > s) new_ed_clusters++; */
                /* clusters++; */

                /* // write new alignment to VCF */
                /* int ref_idx = 0, alt_idx = 0; */
                /* for(size_t cig_idx = 0; cig_idx < cig.size();) { */
                /*     int indel_len = 0; */
                /*     switch (cig[cig_idx]) { */

                /*         case PTR_DIAG: // no variant, update pointers */
                /*             cig_idx += 2; */
                /*             ref_idx++; */
                /*             alt_idx++; */
                /*             break; */

                /*         case PTR_LEFT: // deletion */
                /*             cig_idx++; indel_len++; */

                /*             // multi-base deletion */
                /*             while (cig_idx < cig.size() && cig[cig_idx] == PTR_LEFT) { */
                /*                 cig_idx++; indel_len++; */
                /*             } */
                /*             results.hapcalls[h][ctg].add_var(beg+ref_idx, */
                /*                     indel_len, h, TYPE_DEL, */ 
                /*                     ref_str.substr(ref_idx, indel_len), */
                /*                     "", 60, 60); */
                /*             ref_idx += indel_len; */
                /*             break; */

                /*         case PTR_UP: // insertion */
                /*             cig_idx++; indel_len++; */

                /*             // multi-base insertion */
                /*             while (cig_idx < cig.size() && cig[cig_idx] == PTR_UP) { */
                /*                 cig_idx++; indel_len++; */
                /*             } */
                /*             results.hapcalls[h][ctg].add_var(beg+ref_idx, */
                /*                     0, h, TYPE_INS, "", */ 
                /*                     alt_str.substr(alt_idx, indel_len), 60, 60); */
                /*             alt_idx += indel_len; */
                /*             break; */
                /*     } */
                /* } */


            } // cluster
        } // contig
    } // hap
    INFO("Edit dist reduced in %i of %i clusters, from %i to %i.", 
            new_ed_clusters, clusters, old_ed, new_ed);

    return results;

}

/******************************************************************************/

clusterData edit_dist(vcfData* calls, vcfData* truth, fastaData* ref) {

    // iterate over each contig
    int distance = 0;
    clusterData clusters(calls->contigs);
    for (std::string ctg : calls->contigs) {

        // initialize variant call info
        variantCalls* cal1_vars = &calls->hapcalls[0][ctg];
        variantCalls* cal2_vars = &calls->hapcalls[1][ctg];
        variantCalls* hap1_vars;
        variantCalls* hap2_vars;
        try {
            hap1_vars = &truth->hapcalls[0][ctg];
            hap2_vars = &truth->hapcalls[1][ctg];
        } catch (const std::exception & e) {
            ERROR("truth VCF does not contain contig '%s'", ctg.data());
        }

        // for each cluster of variants (merge calls and truth haps)
        size_t cal1_clust_beg_idx = 0;
        size_t cal2_clust_beg_idx = 0;
        size_t hap1_clust_beg_idx = 0;
        size_t hap2_clust_beg_idx = 0;
        while (cal1_clust_beg_idx < cal1_vars->clusters.size()-1 || // last cluster added for end
               cal2_clust_beg_idx < cal2_vars->clusters.size()-1 ||
               hap1_clust_beg_idx < hap1_vars->clusters.size()-1 ||
               hap2_clust_beg_idx < hap2_vars->clusters.size()-1) {

            // set all start positions
            size_t cal1_clust_end_idx = cal1_clust_beg_idx;
            size_t cal2_clust_end_idx = cal2_clust_beg_idx;
            size_t hap1_clust_end_idx = hap1_clust_beg_idx;
            size_t hap2_clust_end_idx = hap2_clust_beg_idx;

            int cal1_pos = (cal1_clust_beg_idx < cal1_vars->clusters.size()-1) ? 
                    cal1_vars->poss[cal1_vars->clusters[cal1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int cal2_pos = (cal2_clust_beg_idx < cal2_vars->clusters.size()-1) ? 
                    cal2_vars->poss[cal2_vars->clusters[cal2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int hap1_pos = (hap1_clust_beg_idx < hap1_vars->clusters.size()-1) ? 
                    hap1_vars->poss[hap1_vars->clusters[hap1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            int hap2_pos = (hap2_clust_beg_idx < hap2_vars->clusters.size()-1) ? 
                    hap2_vars->poss[hap2_vars->clusters[hap2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();

            // initialize cluster merging with first to start
            int curr_end_pos = 0;
            if (cal1_pos <= hap1_pos && cal1_pos <= cal2_pos && cal1_pos <= hap2_pos) {
                cal1_clust_end_idx += 1;
                curr_end_pos = cal1_vars->poss[cal1_vars->clusters[cal1_clust_end_idx]-1] +
                        cal1_vars->rlens[cal1_vars->clusters[cal1_clust_end_idx]-1] + 1;
                cal1_pos = (cal1_clust_end_idx < cal1_vars->clusters.size()-1) ? 
                    cal1_vars->poss[cal1_vars->clusters[cal1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (cal2_pos <= hap1_pos && cal2_pos <= cal1_pos && 
                    cal2_pos <= hap2_pos) {
                cal2_clust_end_idx += 1;
                curr_end_pos = cal2_vars->poss[cal2_vars->clusters[cal2_clust_end_idx]-1] +
                        cal2_vars->rlens[cal2_vars->clusters[cal2_clust_end_idx]-1] + 1;
                cal2_pos = (cal2_clust_end_idx < cal2_vars->clusters.size()-1) ? 
                    cal2_vars->poss[cal2_vars->clusters[cal2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else if (hap1_pos <= cal1_pos && hap1_pos <= cal2_pos && 
                    hap1_pos <= hap2_pos) {
                hap1_clust_end_idx += 1;
                curr_end_pos = hap1_vars->poss[hap1_vars->clusters[hap1_clust_end_idx]-1] +
                        hap1_vars->rlens[hap1_vars->clusters[hap1_clust_end_idx]-1] + 1;
                hap1_pos = (hap1_clust_end_idx < hap1_vars->clusters.size()-1) ? 
                    hap1_vars->poss[hap1_vars->clusters[hap1_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            } 
            else {
                hap2_clust_end_idx += 1;
                curr_end_pos = hap2_vars->poss[hap2_vars->clusters[hap2_clust_end_idx]-1] +
                        hap2_vars->rlens[hap2_vars->clusters[hap2_clust_end_idx]-1] + 1;
                hap2_pos = (hap2_clust_end_idx < hap2_vars->clusters.size()-1) ? 
                    hap2_vars->poss[hap2_vars->clusters[hap2_clust_end_idx]]-1 : 
                    std::numeric_limits<int>::max();
            }


            // keep expanding cluster while possible
            bool just_merged = true;
            while (just_merged) {
                just_merged = false;
                while (hap1_pos < curr_end_pos + g.gap) {
                    hap1_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            hap1_vars->poss[hap1_vars->clusters[hap1_clust_end_idx]-1] + 
                            hap1_vars->rlens[hap1_vars->clusters[hap1_clust_end_idx]-1] + 1);
                    hap1_pos = (hap1_clust_end_idx < hap1_vars->clusters.size()-1) ? 
                        hap1_vars->poss[hap1_vars->clusters[hap1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (hap2_pos < curr_end_pos + g.gap) {
                    hap2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            hap2_vars->poss[hap2_vars->clusters[hap2_clust_end_idx]-1] + 
                            hap2_vars->rlens[hap2_vars->clusters[hap2_clust_end_idx]-1] + 1);
                    hap2_pos = (hap2_clust_end_idx < hap2_vars->clusters.size()-1) ? 
                        hap2_vars->poss[hap2_vars->clusters[hap2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (cal1_pos < curr_end_pos + g.gap) {
                    cal1_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos, 
                            cal1_vars->poss[cal1_vars->clusters[cal1_clust_end_idx]-1] + 
                            cal1_vars->rlens[cal1_vars->clusters[cal1_clust_end_idx]-1] + 1);
                    cal1_pos = (cal1_clust_end_idx < cal1_vars->clusters.size()-1) ? 
                        cal1_vars->poss[cal1_vars->clusters[cal1_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
                while (cal2_pos < curr_end_pos + g.gap) {
                    cal2_clust_end_idx += 1;
                    curr_end_pos = std::max(curr_end_pos,
                            cal2_vars->poss[cal2_vars->clusters[cal2_clust_end_idx]-1] + 
                            cal2_vars->rlens[cal2_vars->clusters[cal2_clust_end_idx]-1] + 1);
                    cal2_pos = (cal2_clust_end_idx < cal2_vars->clusters.size()-1) ? 
                        cal2_vars->poss[cal2_vars->clusters[cal2_clust_end_idx]]-1 : 
                        std::numeric_limits<int>::max();
                    just_merged = true;
                }
            }

            // get supercluster start/end positions (allowing empty haps)
            int beg_pos = std::numeric_limits<int>::max();
            int end_pos = -1;
            if (cal1_clust_end_idx - cal1_clust_beg_idx) { // cal1 vars present
                beg_pos = std::min(beg_pos, 
                        cal1_vars->poss[cal1_vars->clusters[cal1_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        cal1_vars->poss[cal1_vars->clusters[cal1_clust_end_idx]-1] + 
                        cal1_vars->rlens[cal1_vars->clusters[cal1_clust_end_idx]-1]+1);
            }
            if (cal2_clust_end_idx - cal2_clust_beg_idx) { // cal2 vars present
                beg_pos = std::min(beg_pos, 
                        cal2_vars->poss[cal2_vars->clusters[cal2_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        cal2_vars->poss[cal2_vars->clusters[cal2_clust_end_idx]-1] + 
                        cal2_vars->rlens[cal2_vars->clusters[cal2_clust_end_idx]-1]+1);
            }
            if (hap1_clust_end_idx - hap1_clust_beg_idx) { // hap1 vars present
                beg_pos = std::min(beg_pos, 
                        hap1_vars->poss[hap1_vars->clusters[hap1_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        hap1_vars->poss[hap1_vars->clusters[hap1_clust_end_idx]-1] + 
                        hap1_vars->rlens[hap1_vars->clusters[hap1_clust_end_idx]-1]+1);
            }
            if (hap2_clust_end_idx - hap2_clust_beg_idx) { // hap2 vars present
                beg_pos = std::min(beg_pos, 
                        hap2_vars->poss[hap2_vars->clusters[hap2_clust_beg_idx]]-1);
                end_pos = std::max(end_pos, 
                        hap2_vars->poss[hap2_vars->clusters[hap2_clust_end_idx]-1] + 
                        hap2_vars->rlens[hap2_vars->clusters[hap2_clust_end_idx]-1]+1);
            }

            // generate ref string
            std::string ref_str = ref->fasta.at(ctg).substr(beg_pos, end_pos-beg_pos);

            // generate cal1 and cal2 strings and pointers
            int cal1_var_idx = cal1_vars->clusters[cal1_clust_beg_idx];
            int cal2_var_idx = cal2_vars->clusters[cal2_clust_beg_idx];
            std::string cal1 = "", cal2 = "", cal1_str = "", cal2_str = ""; 
            std::vector<int> cal1_ptrs, cal2_ptrs;
            for (int cal1_ref_pos = beg_pos, cal2_ref_pos = beg_pos; cal1_ref_pos < end_pos || cal2_ref_pos < end_pos; ) {

                // CONSIDER HAP1 ONLY, PRIOR REFERENCE POSITION
                if (cal1_ref_pos < cal2_ref_pos) {
                    if (cal1_ref_pos == cal1_vars->poss[cal1_var_idx]) { // in cal1 variant
                        switch (cal1_vars->types[cal1_var_idx]) {
                            case TYPE_INS:
                                cal1 += cal1_vars->alts[cal1_var_idx];
                                cal1_ptrs.insert(cal1_ptrs.end(), cal1_vars->alts[cal1_var_idx].size(), -1);
                                cal1_str += GREEN(cal1_vars->alts[cal1_var_idx]);
                                break;
                            case TYPE_DEL:
                                cal1_str += std::string(cal1_vars->refs[cal1_var_idx].size(), ' ');
                                cal1_ref_pos += cal1_vars->refs[cal1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                cal1 += cal1_vars->alts[cal1_var_idx];
                                cal1_ptrs.insert(cal1_ptrs.end(), cal1_vars->alts[cal1_var_idx].size(), -1);
                                cal1_str += GREEN(cal1_vars->alts[cal1_var_idx]);
                                cal1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                cal1 += cal1_vars->alts[cal1_var_idx];
                                cal1_ptrs.insert(cal1_ptrs.end(), cal1_vars->alts[cal1_var_idx].size(), -1);
                                cal1_str += std::string(cal1_vars->refs[cal1_var_idx].size(), ' ') 
                                    + GREEN(cal1_vars->alts[cal1_var_idx]);
                                cal1_ref_pos += cal1_vars->refs[cal1_var_idx].size();
                                break;
                        }
                        cal1_var_idx++; // next variant
                    } else { // no cal1 variant, in cal2 variant
                        try {
                            cal1 += ref->fasta.at(ctg)[cal1_ref_pos];
                            cal1_ptrs.push_back(-1);
                            cal1_str += ref->fasta.at(ctg)[cal1_ref_pos];
                            cal1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // CONSIDER HAP2 ONLY, PRIOR REFERENCE POSITION
                else if (cal2_ref_pos < cal1_ref_pos) {
                    if (cal2_ref_pos == cal2_vars->poss[cal2_var_idx]) { // in cal2 variant
                        switch (cal2_vars->types[cal2_var_idx]) {
                            case TYPE_INS:
                                cal2 += cal2_vars->alts[cal2_var_idx];
                                cal2_ptrs.insert(cal2_ptrs.end(), cal2_vars->alts[cal2_var_idx].size(), -1);
                                cal2_str += GREEN(cal2_vars->alts[cal2_var_idx]);
                                break;
                            case TYPE_DEL:
                                cal2_str += std::string(cal2_vars->refs[cal2_var_idx].size(), ' ');
                                cal2_ref_pos += cal2_vars->refs[cal2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                cal2 += cal2_vars->alts[cal2_var_idx];
                                cal2_ptrs.insert(cal2_ptrs.end(), cal2_vars->alts[cal2_var_idx].size(), -1);
                                cal2_str += GREEN(cal2_vars->alts[cal2_var_idx]);
                                cal2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                cal2 += cal2_vars->alts[cal2_var_idx];
                                cal2_ptrs.insert(cal2_ptrs.end(), cal2_vars->alts[cal2_var_idx].size(), -1);
                                cal2_str += std::string(cal2_vars->refs[cal2_var_idx].size(), ' ') 
                                    + GREEN(cal2_vars->alts[cal2_var_idx]);
                                cal2_ref_pos += cal2_vars->refs[cal2_var_idx].size();
                                break;
                        }
                        cal2_var_idx++; // next variant
                    } else { // match
                        try {
                            cal2 += ref->fasta.at(ctg)[cal2_ref_pos];
                            cal2_ptrs.push_back(-1);
                            cal2_str += ref->fasta.at(ctg)[cal2_ref_pos];
                            cal2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // REFERENCE POSITIONS MATCH! POTENTIAL TRANSITIONS
                else {
                    bool cal1_var = false;
                    if (cal1_ref_pos == cal1_vars->poss[cal1_var_idx]) { // in cal1 variant
                        cal1_var = true;
                        switch (cal1_vars->types[cal1_var_idx]) {
                            case TYPE_INS:
                                cal1 += cal1_vars->alts[cal1_var_idx];
                                cal1_ptrs.insert(cal1_ptrs.end(), cal1_vars->alts[cal1_var_idx].size(), -1);
                                cal1_str += GREEN(cal1_vars->alts[cal1_var_idx]);
                                break;
                            case TYPE_DEL:
                                cal1_str += std::string(cal1_vars->refs[cal1_var_idx].size(), ' ');
                                cal1_ref_pos += cal1_vars->refs[cal1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                cal1 += cal1_vars->alts[cal1_var_idx];
                                cal1_ptrs.insert(cal1_ptrs.end(), cal1_vars->alts[cal1_var_idx].size(), -1);
                                cal1_str += GREEN(cal1_vars->alts[cal1_var_idx]);
                                cal1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                cal1 += cal1_vars->alts[cal1_var_idx];
                                cal1_ptrs.insert(cal1_ptrs.end(), cal1_vars->alts[cal1_var_idx].size(), -1);
                                cal1_str += std::string(cal1_vars->refs[cal1_var_idx].size(), ' ') 
                                    + GREEN(cal1_vars->alts[cal1_var_idx]);
                                cal1_ref_pos += cal1_vars->refs[cal1_var_idx].size();
                                break;
                        }
                        cal1_var_idx++; // next variant
                    } 

                    bool cal2_var = false;
                    if (cal2_ref_pos == cal2_vars->poss[cal2_var_idx]) { // in cal2 variant
                        cal2_var = true;
                        switch (cal2_vars->types[cal2_var_idx]) {
                            case TYPE_INS:
                                cal2 += cal2_vars->alts[cal2_var_idx];
                                cal2_ptrs.insert(cal2_ptrs.end(), cal2_vars->alts[cal2_var_idx].size(), -1);
                                cal2_str += GREEN(cal2_vars->alts[cal2_var_idx]);
                                break;
                            case TYPE_DEL:
                                cal2_str += std::string(cal2_vars->refs[cal2_var_idx].size(), ' ');
                                cal2_ref_pos += cal2_vars->refs[cal2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                cal2 += cal2_vars->alts[cal2_var_idx];
                                cal2_ptrs.insert(cal2_ptrs.end(), cal2_vars->alts[cal2_var_idx].size(), -1);
                                cal2_str += GREEN(cal2_vars->alts[cal2_var_idx]);
                                cal2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                cal2 += cal2_vars->alts[cal2_var_idx];
                                cal2_ptrs.insert(cal2_ptrs.end(), cal2_vars->alts[cal2_var_idx].size(), -1);
                                cal2_str += std::string(cal2_vars->refs[cal2_var_idx].size(), ' ') 
                                    + GREEN(cal2_vars->alts[cal2_var_idx]);
                                cal2_ref_pos += cal2_vars->refs[cal2_var_idx].size();
                                break;
                        }
                        cal2_var_idx++; // next variant

                    } 
                    
                    // ONE HAPLOTYPE WAS A VARIANT, INVALID POINTERS
                    if (!cal1_var && cal2_var) {
                        try {
                            cal1 += ref->fasta.at(ctg)[cal1_ref_pos];
                            cal1_ptrs.push_back(-1);
                            cal1_str += ref->fasta.at(ctg)[cal1_ref_pos];
                            cal1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                    if (cal1_var && !cal2_var) {
                        try {
                            cal2 += ref->fasta.at(ctg)[cal2_ref_pos];
                            cal2_ptrs.push_back(-1);
                            cal2_str += ref->fasta.at(ctg)[cal2_ref_pos];
                            cal2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }

                    // BOTH MATCH REFERENCE, ADD POINTERS
                    if (!cal1_var && !cal2_var) { // add pointers
                        try {
                            cal1_ptrs.push_back(cal2.size());
                            cal2_ptrs.push_back(cal1.size());
                            cal1 += ref->fasta.at(ctg)[cal1_ref_pos];
                            cal1_str += BLUE(ref->fasta.at(ctg)[cal1_ref_pos]);
                            cal1_ref_pos++;
                            cal2 += ref->fasta.at(ctg)[cal2_ref_pos];
                            cal2_str += BLUE(ref->fasta.at(ctg)[cal2_ref_pos]);
                            cal2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }
            }

            // generate hap1 and hap2 strings and pointers
            int hap1_var_idx = hap1_vars->clusters[hap1_clust_beg_idx];
            int hap2_var_idx = hap2_vars->clusters[hap2_clust_beg_idx];
            std::string hap1 = "", hap2 = "", hap1_str = "", hap2_str = ""; 
            std::vector<int> hap1_ptrs, hap2_ptrs;
            for (int hap1_ref_pos = beg_pos, hap2_ref_pos = beg_pos; hap1_ref_pos < end_pos || hap2_ref_pos < end_pos; ) {

                // CONSIDER HAP1 ONLY, PRIOR REFERENCE POSITION
                if (hap1_ref_pos < hap2_ref_pos) {
                    if (hap1_ref_pos == hap1_vars->poss[hap1_var_idx]) { // in hap1 variant
                        switch (hap1_vars->types[hap1_var_idx]) {
                            case TYPE_INS:
                                hap1 += hap1_vars->alts[hap1_var_idx];
                                hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                                hap1_str += GREEN(hap1_vars->alts[hap1_var_idx]);
                                break;
                            case TYPE_DEL:
                                hap1_str += std::string(hap1_vars->refs[hap1_var_idx].size(), ' ');
                                hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                hap1 += hap1_vars->alts[hap1_var_idx];
                                hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                                hap1_str += GREEN(hap1_vars->alts[hap1_var_idx]);
                                hap1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                hap1 += hap1_vars->alts[hap1_var_idx];
                                hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                                hap1_str += std::string(hap1_vars->refs[hap1_var_idx].size(), ' ') 
                                    + GREEN(hap1_vars->alts[hap1_var_idx]);
                                hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                                break;
                        }
                        hap1_var_idx++; // next variant
                    } else { // no hap1 variant, in hap2 variant
                        try {
                            hap1 += ref->fasta.at(ctg)[hap1_ref_pos];
                            hap1_ptrs.push_back(-1);
                            hap1_str += ref->fasta.at(ctg)[hap1_ref_pos];
                            hap1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // CONSIDER HAP2 ONLY, PRIOR REFERENCE POSITION
                else if (hap2_ref_pos < hap1_ref_pos) {
                    if (hap2_ref_pos == hap2_vars->poss[hap2_var_idx]) { // in hap2 variant
                        switch (hap2_vars->types[hap2_var_idx]) {
                            case TYPE_INS:
                                hap2 += hap2_vars->alts[hap2_var_idx];
                                hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                                hap2_str += GREEN(hap2_vars->alts[hap2_var_idx]);
                                break;
                            case TYPE_DEL:
                                hap2_str += std::string(hap2_vars->refs[hap2_var_idx].size(), ' ');
                                hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                hap2 += hap2_vars->alts[hap2_var_idx];
                                hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                                hap2_str += GREEN(hap2_vars->alts[hap2_var_idx]);
                                hap2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                hap2 += hap2_vars->alts[hap2_var_idx];
                                hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                                hap2_str += std::string(hap2_vars->refs[hap2_var_idx].size(), ' ') 
                                    + GREEN(hap2_vars->alts[hap2_var_idx]);
                                hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                                break;
                        }
                        hap2_var_idx++; // next variant
                    } else { // match
                        try {
                            hap2 += ref->fasta.at(ctg)[hap2_ref_pos];
                            hap2_ptrs.push_back(-1);
                            hap2_str += ref->fasta.at(ctg)[hap2_ref_pos];
                            hap2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // REFERENCE POSITIONS MATCH! POTENTIAL TRANSITIONS
                else {
                    bool hap1_var = false;
                    if (hap1_ref_pos == hap1_vars->poss[hap1_var_idx]) { // in hap1 variant
                        hap1_var = true;
                        switch (hap1_vars->types[hap1_var_idx]) {
                            case TYPE_INS:
                                hap1 += hap1_vars->alts[hap1_var_idx];
                                hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                                hap1_str += GREEN(hap1_vars->alts[hap1_var_idx]);
                                break;
                            case TYPE_DEL:
                                hap1_str += std::string(hap1_vars->refs[hap1_var_idx].size(), ' ');
                                hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                                break;
                            case TYPE_SUB:
                                hap1 += hap1_vars->alts[hap1_var_idx];
                                hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                                hap1_str += GREEN(hap1_vars->alts[hap1_var_idx]);
                                hap1_ref_pos++;
                                break;
                            case TYPE_GRP:
                                hap1 += hap1_vars->alts[hap1_var_idx];
                                hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                                hap1_str += std::string(hap1_vars->refs[hap1_var_idx].size(), ' ') 
                                    + GREEN(hap1_vars->alts[hap1_var_idx]);
                                hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                                break;
                        }
                        hap1_var_idx++; // next variant
                    } 

                    bool hap2_var = false;
                    if (hap2_ref_pos == hap2_vars->poss[hap2_var_idx]) { // in hap2 variant
                        hap2_var = true;
                        switch (hap2_vars->types[hap2_var_idx]) {
                            case TYPE_INS:
                                hap2 += hap2_vars->alts[hap2_var_idx];
                                hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                                hap2_str += GREEN(hap2_vars->alts[hap2_var_idx]);
                                break;
                            case TYPE_DEL:
                                hap2_str += std::string(hap2_vars->refs[hap2_var_idx].size(), ' ');
                                hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                                break;
                            case TYPE_SUB:
                                hap2 += hap2_vars->alts[hap2_var_idx];
                                hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                                hap2_str += GREEN(hap2_vars->alts[hap2_var_idx]);
                                hap2_ref_pos++;
                                break;
                            case TYPE_GRP:
                                hap2 += hap2_vars->alts[hap2_var_idx];
                                hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                                hap2_str += std::string(hap2_vars->refs[hap2_var_idx].size(), ' ') 
                                    + GREEN(hap2_vars->alts[hap2_var_idx]);
                                hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                                break;
                        }
                        hap2_var_idx++; // next variant

                    } 
                    
                    // ONE HAPLOTYPE WAS A VARIANT, INVALID POINTERS
                    if (!hap1_var && hap2_var) {
                        try {
                            hap1 += ref->fasta.at(ctg)[hap1_ref_pos];
                            hap1_ptrs.push_back(-1);
                            hap1_str += ref->fasta.at(ctg)[hap1_ref_pos];
                            hap1_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                    if (hap1_var && !hap2_var) {
                        try {
                            hap2 += ref->fasta.at(ctg)[hap2_ref_pos];
                            hap2_ptrs.push_back(-1);
                            hap2_str += ref->fasta.at(ctg)[hap2_ref_pos];
                            hap2_ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }

                    // BOTH MATCH REFERENCE, ADD POINTERS
                    if (!hap1_var && !hap2_var) { // add pointers
                        try {
                            hap1_ptrs.push_back(hap2.size());
                            hap2_ptrs.push_back(hap1.size());
                            hap1 += ref->fasta.at(ctg)[hap1_ref_pos];
                            hap1_str += BLUE(ref->fasta.at(ctg)[hap1_ref_pos]);
                            hap1_ref_pos++;
                            hap2 += ref->fasta.at(ctg)[hap2_ref_pos];
                            hap2_str += BLUE(ref->fasta.at(ctg)[hap2_ref_pos]);
                            hap2_ref_pos++;
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
            std::vector<std::string> cals {cal1, cal1, cal2, cal2};
            std::vector<std::string> haps {hap1, hap2, hap1, hap2};
            std::vector<int> call_lens = 
                    {int(cal1.size()), int(cal1.size()), int(cal2.size()), int(cal2.size())};
            std::vector<int> hap_lens = 
                    {int(hap1.size()), int(hap2.size()), int(hap1.size()), int(hap2.size())};
            std::vector< std::vector< std::vector<int> > > offs(4), ptrs(4);

            // for each combination of calls and haps
            for(int i = 0; i < 4; i++) {

                int mat_len = call_lens[i] + hap_lens[i] - 1;
                offs[i].push_back(std::vector<int>(mat_len,-2));
                offs[i][0][call_lens[i]-1] = -1;
                ptrs[i].push_back(std::vector<int>(mat_len,PTR_NONE));
                bool done = false;
                while (true) {

                    // EXTEND WAVEFRONT
                    for (int d = 0; d < mat_len; d++) {
                        int off = offs[i][s[i]][d];
                        int diag = d + 1 - call_lens[i];

                        // don't allow starting from untouched cells
                        if (off == -2) continue;

                        // check that it's within matrix
                        if (diag + off + 1 < 0) continue;
                        if (off > call_lens[i] - 1) continue;
                        if (diag + off > hap_lens[i] - 1) continue;

                        // extend
                        while (off < call_lens[i] - 1 && 
                               diag + off < hap_lens[i] - 1) {
                            if (cals[i][off+1] == haps[i][diag+off+1]) off++;
                            else break;
                        }
                        offs[i][s[i]][d] = off;

                        // finish if done
                        if (off == call_lens[i] - 1 && 
                            off + diag == hap_lens[i] - 1)
                        { done = true; break; }

                    }
                    if (done) break;


                    // NEXT WAVEFRONT
                    // add wavefront, fill edge cells
                    offs[i].push_back(std::vector<int>(mat_len, -2));
                    // bottom left cells
                    if (s[i]+1 == call_lens[i]-1)
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
                        int offtop = (offs[i][s[i]][d+1] == -2) ? 
                            -2 : offs[i][s[i]][d+1]+1;
                        if (offleft >= offtop) {
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
            int dist = std::min(s[CAL1_HAP1]+s[CAL2_HAP2], 
                    s[CAL2_HAP1]+s[CAL1_HAP2]);
            if (g.print_verbosity >= 1 && dist) {
                // print cluster info
                printf("\n\nCAL1: %zu clusters\n", cal1_clust_end_idx-cal1_clust_beg_idx);
                for(size_t i = cal1_clust_beg_idx; i < cal1_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, cal1_vars->clusters[i+1]-cal1_vars->clusters[i]);
                    for(int j = cal1_vars->clusters[i]; j < cal1_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), cal1_vars->poss[j], 
                                cal1_vars->refs[j].size() ? cal1_vars->refs[j].data() : "_", 
                                cal1_vars->alts[j].size() ? cal1_vars->alts[j].data() : "_");
                    }
                }
                printf("CAL2: %zu clusters\n", cal2_clust_end_idx-cal2_clust_beg_idx);
                for(size_t i = cal2_clust_beg_idx; i < cal2_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, cal2_vars->clusters[i+1]-cal2_vars->clusters[i]);
                    for(int j = cal2_vars->clusters[i]; j < cal2_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), cal2_vars->poss[j], 
                                cal2_vars->refs[j].size() ? cal2_vars->refs[j].data() : "_", 
                                cal2_vars->alts[j].size() ? cal2_vars->alts[j].data() : "_");
                    }
                }
                printf("HAP1: %zu clusters\n", hap1_clust_end_idx-hap1_clust_beg_idx);
                for(size_t i = hap1_clust_beg_idx; i < hap1_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, hap1_vars->clusters[i+1]-hap1_vars->clusters[i]);
                    for(int j = hap1_vars->clusters[i]; j < hap1_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), hap1_vars->poss[j], 
                                hap1_vars->refs[j].size() ? hap1_vars->refs[j].data() : "_", 
                                hap1_vars->alts[j].size() ? hap1_vars->alts[j].data() : "_");
                    }
                }
                printf("HAP2: %zu clusters\n", hap2_clust_end_idx-hap2_clust_beg_idx);
                for(size_t i = hap2_clust_beg_idx; i < hap2_clust_end_idx; i++) {
                    printf("\tGroup %zu: %d variants\n", i, hap2_vars->clusters[i+1]-hap2_vars->clusters[i]);
                    for(int j = hap2_vars->clusters[i]; j < hap2_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), hap2_vars->poss[j], 
                                hap2_vars->refs[j].size() ? hap2_vars->refs[j].data() : "_", 
                                hap2_vars->alts[j].size() ? hap2_vars->alts[j].data() : "_");
                    }
                }
                
                printf("ORIG: %s\n", ref_str.data());
                printf("CAL1: %s\n", cal1_str.data());
                printf("CAL2: %s\n", cal2_str.data());
                printf("HAP1: %s\n", hap1_str.data());
                printf("HAP2: %s\n", hap2_str.data());
            }

            // DEBUG PRINT
            if (g.print_verbosity >= 2 && dist) {
                for(int h = 0; h < 4; h++) { // 4 alignments
                    printf("\n%s ALIGNMENT (distance %d)\n", 
                            aln_strs[h].data(), s[h]);

                    // create array
                    std::vector< std::vector<char> > ptr_str;
                    for (int i = 0; i < call_lens[h]; i++)
                        ptr_str.push_back(std::vector<char>(hap_lens[h], '.'));

                    // modify array with pointers
                    int mat_len = call_lens[h] + hap_lens[h] - 1;
                    for (int si = 0; si <= s[h]; si++) {
                        for(int di = 0; di < mat_len; di++) {
                            int diag = di + 1 - call_lens[h];
                            int off = offs[h][si][di];

                            // check that indices are within bounds
                            int cal_pos = off;
                            int hap_pos = diag + off;
                            if (cal_pos < 0 || hap_pos < 0) continue;
                            if (cal_pos > call_lens[h]-1 || 
                                    hap_pos > hap_lens[h]-1) continue;

                            // special case: main diag, no previous edits
                            if (si == 0 && diag == 0) {
                                while (cal_pos >= 0) {
                                    ptr_str[cal_pos--][hap_pos--] = '\\';
                                }
                            }
                            // left edge
                            else if (cal_pos == 0) {
                                ptr_str[cal_pos][hap_pos] = '-';
                            } 
                            // top edge
                            else if (hap_pos == 0) {
                                ptr_str[cal_pos][hap_pos] = '|';
                            } 
                            else {
                                // follow diagonal
                                int top_off = (di < mat_len-1) ? offs[h][si-1][di+1]+1 : -2;
                                int left_off = (di > 0) ? offs[h][si-1][di-1] : -2;
                                while (cal_pos > 0 && hap_pos > 0 && 
                                        cal_pos > top_off && cal_pos > left_off) {
                                    ptr_str[cal_pos--][hap_pos--] = '\\';
                                }
                                // check left/up
                                if (cal_pos == top_off) {
                                    ptr_str[cal_pos][hap_pos] = '|';
                                } else if (cal_pos == left_off) {
                                    ptr_str[cal_pos][hap_pos] = '-';
                                }
                            }
                        }
                    }

                    // print array
                    for (int i = -1; i < call_lens[h]; i++) {
                        for (int j = -1; j < hap_lens[h]; j++) {
                            if (i < 0 && j < 0) {
                                printf("  ");
                            }
                            else if (i < 0) {
                                printf("%c", haps[h][j]);
                            } else if (j < 0) {
                                printf("\n%c ", cals[h][i]);
                            } else {
                                printf("%c", ptr_str[i][j]);
                            }
                        }
                    }
                    printf("\n");

                } // 4 alignments
            } // debug print

            // get cluster phasing
            int orig_phase_dist = s[CAL1_HAP1] + s[CAL2_HAP2];
            int swap_phase_dist = s[CAL2_HAP1] + s[CAL1_HAP2];
            int phase = PHASE_NONE; // default either way
            if (orig_phase_dist < swap_phase_dist) phase = PHASE_ORIG;
            if (swap_phase_dist < orig_phase_dist) phase = PHASE_SWAP;

            // save alignment information
            clusters.clusters[ctg].add(
                    cal1_vars, cal1_clust_beg_idx, cal1_clust_end_idx,
                    cal2_vars, cal2_clust_beg_idx, cal2_clust_end_idx,
                    hap1_vars, hap1_clust_beg_idx, hap1_clust_end_idx,
                    hap2_vars, hap2_clust_beg_idx, hap2_clust_end_idx,
                    phase, orig_phase_dist, swap_phase_dist);

            // reset for next merged cluster
            cal1_clust_beg_idx = cal1_clust_end_idx;
            cal2_clust_beg_idx = cal2_clust_end_idx;
            hap1_clust_beg_idx = hap1_clust_end_idx;
            hap2_clust_beg_idx = hap2_clust_end_idx;
            distance += std::min(orig_phase_dist, swap_phase_dist);

        } // each cluster
    } // each contig
    INFO("Total edit distance: %d", distance);

    return clusters;
}
