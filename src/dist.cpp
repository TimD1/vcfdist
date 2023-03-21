#include <string>
#include <set>
#include <vector>
#include <cstdio>
#include <chrono>
#include <utility>
#include <queue>

#include "dist.h"
#include "edit.h"
#include "print.h"
#include "cluster.h"

template <typename T>
bool contains(const std::unordered_set<T> & wave, const T & idx) {
    return wave.find(idx) != wave.end();
}

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
                        default: ERROR("Unexpected variant type (%i)", vars->types[var]);
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
                            default: 
                                ERROR("Unexpected variant type (%i)", vars->types[var]);
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
                            ERROR("Contig '%s' not present in reference FASTA (edit_dist_realign)",
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
                std::queue<idx2> curr_q, next_q;

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
                        idx2 cell = curr_q.front(); curr_q.pop();
                        int alt_idx = cell.qi;
                        int ref_idx = cell.ri;

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
                        idx2 cell = next_q.front(); next_q.pop();
                        int alt_idx = cell.qi;
                        int ref_idx = cell.ri;

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
                    ptrs[alt_idx][ref_idx] |= PTR_LPATH; // color print path
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
                if (g.verbosity > 1 && old_edits_cluster != new_edits_cluster) {
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
                                    TYPE_SUB, BED_INSIDE, std::string(1,ref_str[ref_idx]), 
                                    std::string(1,alt_str[alt_idx]), 
                                    GT_REF_REF, g.max_qual, g.max_qual);
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
                                    indel_len, h, TYPE_DEL, BED_INSIDE,
                                    ref_str.substr(ref_idx, indel_len),
                                    "", GT_REF_REF, g.max_qual, g.max_qual);
                            ref_idx += indel_len;
                            break;

                        case PTR_UP: // insertion
                            cig_idx++; indel_len++;

                            // multi-base insertion
                            while (cig_idx < cig.size() && cig[cig_idx] == PTR_UP) {
                                cig_idx++; indel_len++;
                            }
                            results.ctg_variants[h][ctg]->add_var(beg+ref_idx,
                                    0, h, TYPE_INS, BED_INSIDE, "", 
                                    alt_str.substr(alt_idx, indel_len), 
                                    GT_REF_REF, g.max_qual, g.max_qual);
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


/* Reverse two strings and their associated pointers for reverse alignment. */
void reverse_ptrs_strs(
        std::string & query, std::string & ref,
        std::vector< std::vector<int> > & query_ptrs, 
        std::vector< std::vector<int> > & ref_ptrs) {

    // reverse strings and vectors
    std::reverse(query.begin(), query.end());
    std::reverse(ref.begin(), ref.end());
    std::reverse(query_ptrs[PTRS].begin(), query_ptrs[PTRS].end());
    std::reverse(query_ptrs[FLAGS].begin(), query_ptrs[FLAGS].end());
    std::reverse(ref_ptrs[PTRS].begin(), ref_ptrs[PTRS].end());
    std::reverse(ref_ptrs[FLAGS].begin(), ref_ptrs[FLAGS].end());

    // update pointers
    for (int & ptr : query_ptrs[PTRS]) { ptr = ref.size()-1 - ptr; }
    for (int & ptr : ref_ptrs[PTRS]) { ptr = query.size()-1 - ptr; }
    for (int & flags : query_ptrs[FLAGS]) { // flip beg/end if only one is set
        if (flags & PTR_VAR_BEG && !(flags & PTR_VAR_END)) {
            flags = PTR_VARIANT | PTR_VAR_END;
        } else if (flags & PTR_VAR_END && !(flags & PTR_VAR_BEG)) {
            flags = PTR_VARIANT | PTR_VAR_BEG;
        }
    }
    for (int & flags : ref_ptrs[FLAGS]) { // flip beg/end if only one is set
        if (flags & PTR_VAR_BEG && !(flags & PTR_VAR_END)) {
            flags = PTR_VARIANT | PTR_VAR_END;
        } else if (flags & PTR_VAR_END && !(flags & PTR_VAR_BEG)) {
            flags = PTR_VARIANT | PTR_VAR_BEG;
        }
    }
    return;
}


/******************************************************************************/


/* Generate the new sequence by applying variants to the reference. */
std::string generate_str(
        std::shared_ptr<fastaData> ref, 
        std::shared_ptr<ctgVariants> vars, const std::string & ctg,
        int var_idx, int end_idx, int beg_pos, int end_pos, int min_qual=0) {

    std::string str = "";
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
                str += ref->fasta.at(ctg).substr(ref_pos, ref_end-ref_pos);
                ref_pos = ref_end;
            } catch (const std::out_of_range & e) {
                ERROR("Contig %s not present in reference FASTA (generate_str)",
                        ctg.data());
                exit(1);
            }
        }
    }
    return str;
}


/******************************************************************************/


/* For each position on hap1 and hap2, it generates pointers from the query to 
 * the reference and vice versa, as well as the strings for alignment.
 */
void generate_ptrs_strs(
        std::string & query_str, std::string & ref_str,
        std::vector< std::vector<int> > & query_ptrs, 
        std::vector< std::vector<int> > & ref_ptrs,
        std::shared_ptr<ctgVariants> query_vars,
        size_t query_clust_beg_idx, size_t query_clust_end_idx,
        int beg_pos, int end_pos, 
        std::shared_ptr<fastaData> ref, const std::string & ctg
        ) {

    // generate query and ref strings and pointers
    query_ptrs.resize(PTR_DIMS);
    ref_ptrs.resize(PTR_DIMS);
    int query_var_idx = query_vars->clusters.size() ? 
            query_vars->clusters[query_clust_beg_idx] : 0;
    int query_end_idx = query_vars->clusters.size() ? 
            query_vars->clusters[query_clust_end_idx] : 0;
    for (int ref_pos = beg_pos; ref_pos < end_pos; ) {

        if (query_var_idx < query_end_idx && 
                ref_pos == query_vars->poss[query_var_idx]) { // start query variant
            switch (query_vars->types[query_var_idx]) {
                case TYPE_INS:
                    query_ptrs[PTRS].insert(query_ptrs[PTRS].end(), 
                            query_vars->alts[query_var_idx].size(), ref_str.size()-1);
                    query_ptrs[FLAGS].insert(query_ptrs[FLAGS].end(),
                            query_vars->alts[query_var_idx].size(), PTR_VARIANT);
                    query_ptrs[FLAGS][ query_ptrs[FLAGS].size()-1] |= PTR_VAR_END;
                    query_ptrs[FLAGS][ query_ptrs[FLAGS].size() - 
                        query_vars->alts[query_var_idx].size()] |= PTR_VAR_BEG;
                    query_str += query_vars->alts[query_var_idx];
                    break;
                case TYPE_DEL:
                    ref_ptrs[PTRS].insert(ref_ptrs[PTRS].end(),
                            query_vars->refs[query_var_idx].size(), query_str.size()-1);
                    ref_ptrs[FLAGS].insert(ref_ptrs[FLAGS].end(),
                            query_vars->refs[query_var_idx].size(), PTR_VARIANT);
                    ref_ptrs[FLAGS][ ref_ptrs[FLAGS].size()-1] |= PTR_VAR_END;
                    ref_ptrs[FLAGS][ ref_ptrs[FLAGS].size() - 
                        query_vars->refs[query_var_idx].size()] |= PTR_VAR_BEG;
                    ref_str += query_vars->refs[query_var_idx];
                    ref_pos += query_vars->refs[query_var_idx].size();
                    break;
                case TYPE_SUB:
                    ref_ptrs[PTRS].push_back(query_str.size());
                    ref_ptrs[FLAGS].push_back(PTR_VARIANT|PTR_VAR_BEG|PTR_VAR_END);
                    query_ptrs[PTRS].push_back(ref_str.size());
                    query_ptrs[FLAGS].push_back(PTR_VARIANT|PTR_VAR_BEG|PTR_VAR_END);
                    ref_str += query_vars->refs[query_var_idx];
                    query_str += query_vars->alts[query_var_idx];
                    ref_pos++;
                    break;
                default:
                    ERROR("Unexpected variant type in generate_ptrs_strs()");
            }
            query_var_idx++; // next variant

        } else { // add all matching ref bases

            try {
                // find next position w/o ref match (next var or end)
                int ref_end = (query_var_idx < query_end_idx) ?
                    query_vars->poss[query_var_idx] : end_pos;

                // add pointers, query
                std::vector<int> new_query_ptrs(ref_end - ref_pos, 0);
                query_ptrs[FLAGS].insert(query_ptrs[FLAGS].end(),
                        new_query_ptrs.begin(), new_query_ptrs.end()); // zeros
                for(size_t i = 0; i < new_query_ptrs.size(); i++)
                    new_query_ptrs[i] = ref_str.size() + i;
                query_ptrs[PTRS].insert(query_ptrs[PTRS].end(), 
                        new_query_ptrs.begin(), new_query_ptrs.end());

                // add pointers, ref
                std::vector<int> new_ref_ptrs(ref_end - ref_pos, 0);
                ref_ptrs[FLAGS].insert(ref_ptrs[FLAGS].end(),
                        new_ref_ptrs.begin(), new_ref_ptrs.end()); // zeros
                for(size_t i = 0; i < new_ref_ptrs.size(); i++)
                    new_ref_ptrs[i] = query_str.size() + i;
                ref_ptrs[PTRS].insert(ref_ptrs[PTRS].end(), 
                        new_ref_ptrs.begin(), new_ref_ptrs.end());

                // add sequence, update positions
                std::string matches = ref->fasta.at(ctg).substr(ref_pos, ref_end-ref_pos);
                query_str += matches;
                ref_str += matches;
                ref_pos = ref_end;

            } catch (const std::out_of_range & e) {
                ERROR("Contig '%s' not present in reference FASTA", ctg.data());
            }
        }
    }
}


/******************************************************************************/


void calc_prec_recall_aln(
        const std::string & query1, const std::string & query2,
        const std::string & truth1, const std::string & truth2, 
        const std::string & ref,
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & ptrs,
        std::vector<int> & pr_query_ref_end
        ) {
    
    // set loop variables
    int ref_len = ref.size();
    std::vector<std::string> query {query1, query1, query2, query2};
    std::vector<std::string> truth {truth1, truth2, truth1, truth2};
    std::vector< std::vector< std::vector<int> > > query_ref_ptrs {
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs {
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs {
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector<int> query_lens = 
            {int(query1.size()), int(query1.size()), int(query2.size()), int(query2.size())};
    std::vector<int> truth_lens = 
            {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};

    std::vector< std::vector< std::vector<bool> > > done;

    // for each combination of query and truth
    for (int i = 0; i < 4; i++) {
        int qi = 2*i + QUERY; // query index (ptrs)
        int ri = 2*i + REF;   // ref index   (ptrs)

        // init full pointer/done matrices
        ptrs.push_back(std::vector< std::vector<int> >(query_lens[i], 
                std::vector<int>(truth_lens[i], PTR_NONE)));
        ptrs.push_back(std::vector< std::vector<int> >(ref_len, 
                std::vector<int>(truth_lens[i], PTR_NONE)));
        ptrs[qi][0][0] = PTR_DIAG;

        done.push_back(std::vector< std::vector<bool> >(query_lens[i], 
                std::vector<bool>(truth_lens[i], false)));
        done.push_back(std::vector< std::vector<bool> >(ref_len, 
                std::vector<bool>(truth_lens[i], false)));
        
        // set first wavefront
        std::queue<idx1> queue; // still to be explored in this wave
        queue.push({qi, 0, 0});
        ptrs[qi][0][0] |= PTR_DIAG;
        done[qi][0][0] = true;
        queue.push({ri, 0, 0});
        ptrs[ri][0][0] |= PTR_DIAG;
        done[ri][0][0] = true;

        // continue looping until full alignment found
        std::unordered_set<idx1> curr_wave; // everything explored this wave
        std::unordered_set<idx1> prev_wave; // everything explored prev wave
        while (true) {
            /* printf("s = %d\n", s[i]); */
            if (queue.empty()) ERROR("Empty queue in 'prec_recall_aln()'.");

            // EXTEND WAVEFRONT (stay at same score)
            while (!queue.empty()) {
                idx1 x = queue.front(); queue.pop();
                /* printf("  x = (%d, %d, %d)\n", x.hi, x.qri, x.ti); */
                prev_wave.insert(x);
                if (x.hi == qi) { // == QUERY
                    // allow match
                    if (x.qri+1 < query_lens[i] && x.ti+1 < truth_lens[i] &&
                            query[i][x.qri+1] == truth[i][x.ti+1]) {
                        if (!done[x.hi][x.qri+1][x.ti+1] && 
                                !contains(curr_wave, idx1(x.hi, x.qri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.qri+1, x.ti+1));
                            curr_wave.insert(idx1(x.hi, x.qri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.qri+1][x.ti+1]) {
                            ptrs[x.hi][x.qri+1][x.ti+1] |= PTR_DIAG;
                        }
                    }
                    // allow phase swaps, outside query/truth variants
                    if ( !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT) && 
                         !(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT)) {
                        if (!done[ri][ query_ref_ptrs[i][PTRS][x.qri] ][x.ti] &&
                                !contains(curr_wave, idx1(ri, query_ref_ptrs[i][PTRS][x.qri], x.ti))) {
                            queue.push(idx1(ri, 
                                    query_ref_ptrs[i][PTRS][x.qri], 
                                    x.ti));
                            curr_wave.insert(idx1(ri, query_ref_ptrs[i][PTRS][x.qri], x.ti));
                        }
                        if (!done[ri][ query_ref_ptrs[i][PTRS][x.qri] ][x.ti]) {
                            ptrs[ri][ query_ref_ptrs[i][PTRS][x.qri] ][x.ti] |= PTR_SWAP;
                        }
                    }
                } else { // x.hi == ri == REF
                    // allow match
                    if (x.qri+1 < ref_len && x.ti+1  < truth_lens[i] &&
                            ref[x.qri+1] == truth[i][x.ti+1]) {
                        if (!done[x.hi][x.qri+1][x.ti+1] &&
                                !contains(curr_wave, idx1(x.hi, x.qri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.qri+1, x.ti+1));
                            curr_wave.insert(idx1(x.hi, x.qri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.qri+1][x.ti+1]) {
                            ptrs[x.hi][x.qri+1][x.ti+1] |= PTR_DIAG;
                        }
                    }
                    // allow phase swaps, outside query/truth variants
                    if ( !(ref_query_ptrs[i][FLAGS][x.qri] & PTR_VARIANT) && 
                         !(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT)) {
                        if (!done[qi][ref_query_ptrs[i][PTRS][x.qri]][x.ti] &&
                                !contains(curr_wave, idx1(qi, ref_query_ptrs[i][PTRS][x.qri], x.ti))) {
                            queue.push(idx1(qi, ref_query_ptrs[i][PTRS][x.qri], x.ti));
                            curr_wave.insert(idx1(qi, ref_query_ptrs[i][PTRS][x.qri], x.ti));
                        }
                        if (!done[qi][ ref_query_ptrs[i][PTRS][x.qri] ][x.ti]) {
                            ptrs[qi][ ref_query_ptrs[i][PTRS][x.qri] ][x.ti] |= PTR_SWAP;
                        }
                    }
                }
            }

            // mark all cells visited this wave as done
            for (auto x : curr_wave) { done[x.hi][x.qri][x.ti] = true; }
            curr_wave.clear();

            // exit if we're done aligning
            if (done[qi][query_lens[i]-1][truth_lens[i]-1] ||
                done[ri][ref_len-1][truth_lens[i]-1]) break;


            // NEXT WAVEFRONT (increase score by one)
            for (auto x : prev_wave) {
                if (x.hi == qi) { // QUERY
                    if (x.qri+1 < query_lens[i]) { // INS
                        if (!done[x.hi][x.qri+1][x.ti] && 
                                !contains(curr_wave,
                                    idx1(x.hi, x.qri+1, x.ti))) {
                            queue.push(idx1(x.hi, x.qri+1, x.ti));
                            curr_wave.insert(idx1(x.hi, x.qri+1, x.ti));
                        }
                        if (!done[x.hi][x.qri+1][x.ti])
                            ptrs[x.hi][x.qri+1][x.ti] |= PTR_UP;
                    }
                    if (x.ti+1 < truth_lens[i]) { // DEL
                        if (!done[x.hi][x.qri][x.ti+1] &&
                                !contains(curr_wave,
                                    idx1(x.hi, x.qri, x.ti+1))) {
                            queue.push(idx1(x.hi, x.qri, x.ti+1));
                            curr_wave.insert(idx1(x.hi, x.qri, x.ti+1));
                        }
                        if (!done[x.hi][x.qri][x.ti+1])
                            ptrs[x.hi][x.qri][x.ti+1] |= PTR_LEFT;
                    }
                    if (x.qri+1 < query_lens[i] && x.ti+1 < truth_lens[i]) { // SUB
                        if (!done[x.hi][x.qri+1][x.ti+1] &&
                                !contains(curr_wave,
                                    idx1(x.hi, x.qri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.qri+1, x.ti+1));
                            curr_wave.insert(idx1(x.hi, x.qri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.qri+1][x.ti+1])
                            ptrs[x.hi][x.qri+1][x.ti+1] |= PTR_SUB;
                    }
                } else { //x.hi == REF
                    if (x.qri+1 < ref_len) { // INS
                        if (!done[x.hi][x.qri+1][x.ti] &&
                                !contains(curr_wave,
                                    idx1(x.hi, x.qri+1, x.ti))) {
                            queue.push(idx1(x.hi, x.qri+1, x.ti));
                            curr_wave.insert(idx1(x.hi, x.qri+1, x.ti));
                        }
                        if (!done[x.hi][x.qri+1][x.ti])
                            ptrs[x.hi][x.qri+1][x.ti] |= PTR_UP;
                    }
                    if (x.ti+1 < truth_lens[i]) { // DEL
                        if (!done[x.hi][x.qri][x.ti+1] &&
                                !contains(curr_wave,
                                    idx1(x.hi, x.qri, x.ti+1))) {
                            queue.push(idx1(x.hi, x.qri, x.ti+1));
                            curr_wave.insert(idx1(x.hi, x.qri, x.ti+1));
                        }
                        if (!done[x.hi][x.qri][x.ti+1])
                            ptrs[x.hi][x.qri][x.ti+1] |= PTR_LEFT;
                    }
                    if (x.qri+1 < ref_len && x.ti+1 < truth_lens[i]) { // SUB
                        if (!done[x.hi][x.qri+1][x.ti+1] &&
                                !contains(curr_wave, 
                                    idx1(x.hi, x.qri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.qri+1, x.ti+1));
                            curr_wave.insert(idx1(x.hi, x.qri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.qri+1][x.ti+1])
                            ptrs[x.hi][x.qri+1][x.ti+1] |= PTR_SUB;
                    }
                }
            }
            prev_wave.clear();
            s[i]++;
        } // while loop (this alignment)

        // save where to start backtrack (prefer ref: omit vars which don't reduce ED)
        if (done[ri][ref_len-1][truth_lens[i]-1]) {
            pr_query_ref_end[i] = ri;
        } else if (done[qi][query_lens[i]-1][truth_lens[i]-1]) {
            pr_query_ref_end[i] = qi;
        } else { ERROR("Alignment not finished in 'prec_recall_aln()'."); }

    } // 4 alignments
} // function


/******************************************************************************/


int store_phase( 
        std::shared_ptr<superclusterData> clusterdata_ptr, 
        const std::string & ctg,
        const std::vector<int> & s
        ) {

    // calculate best phasing
    int orig_phase_dist = s[QUERY1_TRUTH1] + s[QUERY2_TRUTH2];
    int swap_phase_dist = s[QUERY2_TRUTH1] + s[QUERY1_TRUTH2];
    int phase = PHASE_NONE; // default either way if equal dist
    if (orig_phase_dist < swap_phase_dist) phase = PHASE_ORIG;
    if (swap_phase_dist < orig_phase_dist) phase = PHASE_SWAP;

    // save alignment information
    clusterdata_ptr->ctg_superclusters[ctg]->add_phasing(
            phase, orig_phase_dist, swap_phase_dist);
    return phase;
}


/******************************************************************************/


void calc_prec_recall_path(
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, 
        const std::string & ctg,
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_scores, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_end
        ) {

    std::vector< std::vector< std::vector<int> > > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector<int> pr_query_ref_beg(4);
    std::vector< std::vector<bool> > ref_loc_sync;

    for (int i = 0; i < 4; i++) {

        // init
        int qi = i*2 + QUERY;
        int ri = i*2 + REF;
        path.push_back(std::vector<idx1>());
        sync.push_back(std::vector<bool>());
        edits.push_back(std::vector<bool>());
        ref_loc_sync.push_back(std::vector<bool>(int(aln_ptrs[ri].size()), true));
        path_ptrs.push_back(std::vector< std::vector<int> >(aln_ptrs[qi].size(), 
            std::vector<int>(aln_ptrs[qi][0].size(), PTR_NONE)));
        path_ptrs.push_back(std::vector< std::vector<int> >(aln_ptrs[ri].size(), 
            std::vector<int>(aln_ptrs[ri][0].size(), PTR_NONE)));
        path_scores.push_back(std::vector< std::vector<int> >(aln_ptrs[qi].size(), 
            std::vector<int>(aln_ptrs[qi][0].size(), -1)));
        path_scores.push_back(std::vector< std::vector<int> >(aln_ptrs[ri].size(), 
            std::vector<int>(aln_ptrs[ri][0].size(), -1)));

        // backtrack start
        std::queue<idx1> queue;
        int start_hi = pr_query_ref_end[i];
        int start_qri = aln_ptrs[start_hi].size()-1;
        int start_ti = aln_ptrs[start_hi][0].size()-1;
        idx1 start(start_hi, start_qri, start_ti);
        path_ptrs[start_hi][start_qri][start_ti] = PTR_DIAG;
        aln_ptrs[start_hi][start_qri][start_ti] |= MAIN_PATH;
        path_scores[start_hi][start_qri][start_ti] = 0;
        queue.push(start);

        while (!queue.empty()) {
            idx1 x = queue.front(); queue.pop(); // current cell

            // ref locations which consume a query base and aren't 
            // on the main diagonal cannot be sync points
            if (aln_ptrs[x.hi][x.qri][x.ti] & (PTR_DIAG | PTR_SUB | PTR_INS)) {
                if (x.hi == ri && x.qri != truth_ref_ptrs[i][PTRS][x.ti])
                    ref_loc_sync[i][x.qri] = false;
                if (x.hi == qi && !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT) && 
                        query_ref_ptrs[i][PTRS][x.qri] != truth_ref_ptrs[i][PTRS][x.ti])
                    ref_loc_sync[i][query_ref_ptrs[i][PTRS][x.qri]] = false;
            }

            if (aln_ptrs[x.hi][x.qri][x.ti] & PTR_DIAG && // MATCH
                    x.qri > 0 && x.ti > 0) {

                // add to path
                idx1 y = idx1(x.hi, x.qri-1, x.ti-1); // next cell
                aln_ptrs[y.hi][y.qri][y.ti] |= PATH;

                // check for fp
                int is_fp = FALSE;
                if (x.hi == ri && !(ref_query_ptrs[i][FLAGS][y.qri] & PTR_VARIANT) &&
                        (ref_query_ptrs[i][FLAGS][x.qri] & PTR_VARIANT ||
                         ref_query_ptrs[i][PTRS][x.qri] != ref_query_ptrs[i][PTRS][y.qri]+1)) {
                    is_fp = TRUE;
                }

                // update score
                if (path_scores[x.hi][x.qri][x.ti]+is_fp > 
                        path_scores[y.hi][y.qri][y.ti]) {
                    if (y.qri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.qri][y.ti] = PTR_DIAG;
                    path_scores[y.hi][y.qri][y.ti] = 
                        path_scores[x.hi][x.qri][x.ti]+is_fp;
                    queue.push(y);
                }

            }

            if (aln_ptrs[x.hi][x.qri][x.ti] & PTR_SUB && // SUB
                    x.qri > 0 && x.ti > 0) {

                // add to path
                idx1 y = idx1(x.hi, x.qri-1, x.ti-1);
                aln_ptrs[y.hi][y.qri][y.ti] |= PATH;

                // check for fp
                int is_fp = FALSE;
                if (x.hi == ri && !(ref_query_ptrs[i][FLAGS][y.qri] & PTR_VARIANT) &&
                        (ref_query_ptrs[i][FLAGS][x.qri] & PTR_VARIANT ||
                         ref_query_ptrs[i][PTRS][x.qri] != ref_query_ptrs[i][PTRS][y.qri]+1)) {
                    is_fp = TRUE;
                }

                // update score
                if (path_scores[x.hi][x.qri][x.ti]+is_fp > 
                        path_scores[y.hi][y.qri][y.ti]) {
                    if (y.qri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.qri][y.ti] = PTR_SUB;
                    path_scores[y.hi][y.qri][y.ti] = 
                        path_scores[x.hi][x.qri][x.ti]+is_fp;
                    queue.push(y);
                }
            }

            if (aln_ptrs[x.hi][x.qri][x.ti] & PTR_INS && x.qri > 0) { // INS

                // add to path
                idx1 y = idx1(x.hi, x.qri-1, x.ti);
                aln_ptrs[y.hi][y.qri][y.ti] |= PATH;

                // check for fp
                int is_fp = FALSE;
                if (x.hi == ri && !(ref_query_ptrs[i][FLAGS][y.qri] & PTR_VARIANT) &&
                        (ref_query_ptrs[i][FLAGS][x.qri] & PTR_VARIANT ||
                        ref_query_ptrs[i][PTRS][x.qri] != ref_query_ptrs[i][PTRS][y.qri]+1)) {
                    is_fp = TRUE;
                }

                // update score
                if (path_scores[x.hi][x.qri][x.ti]+is_fp > 
                        path_scores[y.hi][y.qri][y.ti]) {
                    if (y.qri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.qri][y.ti] = PTR_INS;
                    path_scores[y.hi][y.qri][y.ti] = 
                        path_scores[x.hi][x.qri][x.ti]+is_fp;
                    queue.push(y);
                } 
            }

            if (aln_ptrs[x.hi][x.qri][x.ti] & PTR_DEL && x.ti > 0) { // DEL

                // add to path
                idx1 y = idx1(x.hi, x.qri, x.ti-1);
                aln_ptrs[y.hi][y.qri][y.ti] |= PATH;

                // update score
                if (path_scores[x.hi][x.qri][x.ti] > 
                        path_scores[y.hi][y.qri][y.ti]) {
                    if (y.qri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.qri][y.ti] = PTR_DEL;
                    path_scores[y.hi][y.qri][y.ti] = 
                        path_scores[x.hi][x.qri][x.ti];
                    queue.push(y);
                }
            }

            if (x.hi == ri) { // REF -> QUERY SWAP
                if (aln_ptrs[x.hi][x.qri][x.ti] & PTR_SWAP &&
                        !(ref_query_ptrs[i][FLAGS][x.qri] & PTR_VARIANT)) {

                    // add to path
                    idx1 y = idx1(qi, ref_query_ptrs[i][PTRS][x.qri], x.ti);
                    aln_ptrs[y.hi][y.qri][y.ti] |= PATH;

                    // update score
                    if (path_scores[x.hi][x.qri][x.ti] > 
                            path_scores[y.hi][y.qri][y.ti]) {
                        if (y.qri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                        path_ptrs  [y.hi][y.qri][y.ti] = PTR_SWAP;
                        path_scores[y.hi][y.qri][y.ti] = 
                            path_scores[x.hi][x.qri][x.ti];
                        queue.push(y);
                    }
                }

            } else { // QUERY -> REF SWAP
                if (aln_ptrs[x.hi][x.qri][x.ti] & PTR_SWAP &&
                        !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT)) {

                    // add to path
                    idx1 y = idx1(ri, query_ref_ptrs[i][PTRS][x.qri], x.ti);
                    aln_ptrs[y.hi][y.qri][y.ti] |= PATH;

                    // update score
                    if (path_scores[x.hi][x.qri][x.ti] > 
                            path_scores[y.hi][y.qri][y.ti]) {
                        if (y.qri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                        path_ptrs  [y.hi][y.qri][y.ti] = PTR_SWAP;
                        path_scores[y.hi][y.qri][y.ti] = 
                            path_scores[x.hi][x.qri][x.ti];
                        queue.push(y);
                    }
                }
            }
        }
    }

    // get path and sync points
    get_prec_recall_path_sync(path, sync, edits, 
            aln_ptrs, path_ptrs, ref_loc_sync,
            query1_ref_ptrs, ref_query1_ptrs, query2_ref_ptrs, ref_query2_ptrs,
            truth1_ref_ptrs, truth2_ref_ptrs, pr_query_ref_beg
    );

}


/******************************************************************************/


void get_prec_recall_path_sync(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        const std::vector< std::vector<bool> > & ref_loc_sync, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_beg
        ) {

    // query <-> ref pointers
    std::vector< std::vector< std::vector<int> > > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };

    for (int i = 0; i < 4; i++) {

        // init
        int qi = i*2 + QUERY;
        int ri = i*2 + REF;

        // path start
        int hi = pr_query_ref_beg[i];
        int qri = 0;
        int ti = 0;

        // first position is sync point
        path[i].push_back(idx1(hi, qri, ti));
        sync[i].push_back(true);
        aln_ptrs[hi][qri][ti] |= MAIN_PATH | PTR_SYNC;
        path_ptrs[hi][qri][ti] |= MAIN_PATH;

        // follow best-path pointers
        while (qri < int(path_ptrs[hi].size()) || ti < int(path_ptrs[hi][0].size())) {
            if (path_ptrs[hi][qri][ti] & PTR_DIAG) {
                qri++; ti++; edits[i].push_back(false);
            } else if (path_ptrs[hi][qri][ti] & PTR_SUB) {
                qri++; ti++; edits[i].push_back(true);
            } else if (path_ptrs[hi][qri][ti] & PTR_INS) {
                qri++; edits[i].push_back(true);
            } else if (path_ptrs[hi][qri][ti] & PTR_DEL) {
                ti++; edits[i].push_back(true);
            } else if (path_ptrs[hi][qri][ti] & PTR_SWAP) {
                if (hi == ri) {
                    hi = qi;
                    qri = ref_query_ptrs[i][PTRS][qri];
                } else if(hi == qi) { // hi == qi
                    hi = ri;
                    qri = query_ref_ptrs[i][PTRS][qri];
                } else {
                    ERROR("Unexpected hap index (%d)", hi);
                } 
                edits[i].push_back(false);
            } else {
                ERROR("No pointer for MAIN_PATH at (%d, %d, %d)", hi, qri, ti);
            }

            if (qri >= int(aln_ptrs[hi].size()) ||
                    ti >= int(aln_ptrs[hi][0].size())) break;

            // add point to path
            path[i].push_back(idx1(hi, qri, ti));
            aln_ptrs[hi][qri][ti] |= MAIN_PATH;
            path_ptrs[hi][qri][ti] |= MAIN_PATH;

            // determine if sync point
            bool is_sync;
            if (hi == qi) {
                if (!(query_ref_ptrs[i][FLAGS][qri] & PTR_VARIANT) &&  // on main diag
                    !(truth_ref_ptrs[i][FLAGS][ti] & PTR_VARIANT) &&
                    query_ref_ptrs[i][PTRS][qri] == truth_ref_ptrs[i][PTRS][ti]) {
                    is_sync = ref_loc_sync[i][ query_ref_ptrs[i][PTRS][qri] ];
                } else {
                    is_sync = false;
                }
            } else { // hi == ri
                if (qri == truth_ref_ptrs[i][PTRS][ti]) { // on main diag
                    is_sync = ref_loc_sync[i][qri];
                } else {
                    is_sync = false;
                }
            }
            if (is_sync) aln_ptrs[hi][qri][ti] |= PTR_SYNC;
            sync[i].push_back(is_sync);
        }
        // last position is sync
        sync[i][sync[i].size()-1] = true;
        idx1 last = path[i][path[i].size()-1];
        aln_ptrs[last.hi][last.qri][last.ti] |= PTR_SYNC;
    }
}


/******************************************************************************/


void calc_prec_recall(
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, 
        const std::string & ctg, const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        const std::vector< std::vector<idx1> > & path,
        const std::vector< std::vector<bool> > & sync,
        const std::vector< std::vector<bool> > & edits,
        const std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        const std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_end, int phase, int print
        ) {

    // set query/truth strings and pointers
    int beg = clusterdata_ptr->ctg_superclusters[ctg]->begs[sc_idx];
    std::vector<std::string> query = {query1, query1, query2, query2};
    std::vector<std::string> truth = {truth1, truth2, truth1, truth2};
    std::vector< std::vector< std::vector<int> > > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };

    // indices into ptr/off matrices depend on decided phasing
    std::vector<int> indices;
    if (phase == PHASE_SWAP) {
        indices.push_back(QUERY1_TRUTH2);
        indices.push_back(QUERY2_TRUTH1);
    } else if (phase == PHASE_ORIG || phase == PHASE_NONE) { // keep
        indices.push_back(QUERY1_TRUTH1);
        indices.push_back(QUERY2_TRUTH2);
    } else {
        ERROR("Unexpected phase (%d)", phase);
    }

    // for only the selected phasing
    for (int i : indices) {

        int ri = i*2 + REF;   // ref index
        int qi = i*2 + QUERY; // call index
        int qhi = i >> 1;      // query hap index
        int thi = i & 1;      // truth hap index

        // set variant ranges
        std::shared_ptr<ctgVariants> query_vars = 
                clusterdata_ptr->ctg_superclusters[ctg]->ctg_variants[QUERY][qhi];
        std::shared_ptr<ctgVariants> truth_vars = 
                clusterdata_ptr->ctg_superclusters[ctg]->ctg_variants[TRUTH][thi];
        int query_beg_idx = query_vars->clusters.size() ? query_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[QUERY][qhi][sc_idx]] : 0;
        int query_end_idx = query_vars->clusters.size() ? query_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[QUERY][qhi][sc_idx+1]] : 0;
        int truth_beg_idx = truth_vars->clusters.size() ? truth_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[TRUTH][thi][sc_idx]] : 0;
        int truth_end_idx = truth_vars->clusters.size() ? truth_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[TRUTH][thi][sc_idx+1]] : 0;

        // debug print
        if (print) printf("Alignment %s, aln_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(aln_ptrs[qi], query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(aln_ptrs[ri], ref, truth[i]);
        if (print) printf("Alignment %s, path_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(path_ptrs[qi], query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(path_ptrs[ri], ref, truth[i]);

        // init
        int hi = pr_query_ref_end[i];
        int qri = aln_ptrs[hi].size()-1;
        int ti = aln_ptrs[hi][0].size()-1;
        int new_ed = 0;
        int query_var_ptr = query_end_idx-1;
        int query_var_pos = query_vars->poss.size() ? query_vars->poss[query_var_ptr] - beg : 0;
        int prev_query_var_ptr = query_var_ptr;
        int truth_var_ptr = truth_end_idx-1;
        int truth_var_pos = truth_vars->poss.size() ? truth_vars->poss[truth_var_ptr] - beg : 0;
        int prev_truth_var_ptr = truth_var_ptr;
        int pidx = path[i].size()-1;

        if (print) printf("\n%s:\n", aln_strs[i].data());
        while (pidx >= 0) {

            // set FP if necessary
            int qref_pos = (hi == ri) ? qri : query_ref_ptrs[i][PTRS][qri];
            if (query_ref_ptrs[i][FLAGS][qri] & PTR_VARIANT) qref_pos = query_var_pos+1;
            while (qref_pos < query_var_pos) {
                if (hi == ri) { // FP
                    query_vars->errtypes[query_var_ptr] = ERRTYPE_FP;
                    query_vars->credit[query_var_ptr] = 0;
                    query_vars->callq[query_var_ptr] = 
                        query_vars->var_quals[query_var_ptr];
                    if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                            query_vars->refs[query_var_ptr].data(),
                            query_vars->alts[query_var_ptr].data(), "FP", 0.0f);
                    query_var_ptr--;
                } else { // passed variant
                    query_var_ptr--;
                }
                query_var_pos = (query_var_ptr < query_beg_idx) ? -1 :
                    query_vars->poss[query_var_ptr] - beg;
            }

            while (!(truth_ref_ptrs[i][FLAGS][ti] & PTR_VARIANT) && 
                    truth_ref_ptrs[i][PTRS][ti] < truth_var_pos) { // passed REF variant
                truth_var_ptr--;
                truth_var_pos = (truth_var_ptr < truth_beg_idx) ? -1 :
                    truth_vars->poss[truth_var_ptr] - beg;
            }

            // sync point: set TP/PP
            if (sync[i][pidx]) {

                // calculate old edit distance
                int old_ed = 0;
                for (int truth_var_idx = prev_truth_var_ptr; 
                        truth_var_idx > truth_var_ptr; truth_var_idx--) {
                    switch (truth_vars->types[truth_var_idx]) {
                        case TYPE_SUB:
                            old_ed += 1;
                            break;
                        case TYPE_INS:
                            old_ed += truth_vars->alts[truth_var_idx].length();
                            break;
                        case TYPE_DEL:
                            old_ed += truth_vars->refs[truth_var_idx].length();
                            break;
                        default:
                            ERROR("Unexpected variant type (%d) in calc_prec_recall().",
                                    truth_vars->types[truth_var_idx]);
                            break;
                    }
                }
                if (old_ed == 0 && truth_var_ptr != prev_truth_var_ptr) 
                    ERROR("Old edit distance 0, variants exist.");

                // get min query var qual in sync section (for truth/query)
                float callq = g.max_qual;
                for (int query_var_idx = prev_query_var_ptr; 
                        query_var_idx > query_var_ptr; query_var_idx--) {
                    callq = std::min(callq, query_vars->var_quals[query_var_idx]);
                }

                // process QUERY variants
                for (int query_var_idx = prev_query_var_ptr; 
                        query_var_idx > query_var_ptr; query_var_idx--) {
                    float credit = 1 - float(new_ed)/old_ed;
                    // don't overwrite FPs
                    if (query_vars->errtypes[query_var_idx] == ERRTYPE_UN) {
                        if (new_ed == 0) { // TP
                            query_vars->errtypes[query_var_idx] = ERRTYPE_TP;
                            query_vars->credit[query_var_idx] = credit;
                            query_vars->callq[query_var_idx] = callq;
                            if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                    query_vars->refs[query_var_idx].data(),
                                    query_vars->alts[query_var_idx].data(), 
                                    "TP", credit);
                        } else if (new_ed == old_ed) { // FP
                            query_vars->errtypes[query_var_idx] = ERRTYPE_FP;
                            query_vars->credit[query_var_idx] = 0;
                            query_vars->callq[query_var_idx] = callq;
                            if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                    query_vars->refs[query_var_idx].data(),
                                    query_vars->alts[query_var_idx].data(), 
                                    "FP", 0.0);
                        } else { // PP
                            query_vars->errtypes[query_var_idx] = ERRTYPE_PP;
                            query_vars->credit[query_var_idx] = credit;
                            query_vars->callq[query_var_idx] = callq;
                            if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                    query_vars->refs[query_var_idx].data(),
                                    query_vars->alts[query_var_idx].data(), 
                                    "PP", credit);
                        }
                    }
                }

                // process TRUTH variants
                for (int truth_var_idx = prev_truth_var_ptr; 
                        truth_var_idx > truth_var_ptr; truth_var_idx--) {
                    float credit = 1 - float(new_ed)/old_ed;
                    if (new_ed == 0) { // TP
                        truth_vars->errtypes[truth_var_idx] = ERRTYPE_TP;
                        truth_vars->credit[truth_var_idx] = credit;
                        truth_vars->callq[truth_var_idx] = callq;
                        if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                truth_vars->refs[truth_var_idx].data(),
                                truth_vars->alts[truth_var_idx].data(), 
                                "TP", credit);
                    } else if (new_ed == old_ed) { // FP call, FN truth
                        truth_vars->errtypes[truth_var_idx] = ERRTYPE_FN;
                        truth_vars->credit[truth_var_idx] = credit;
                        truth_vars->callq[truth_var_idx] = g.max_qual;
                        if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                truth_vars->refs[truth_var_idx].data(),
                                truth_vars->alts[truth_var_idx].data(), 
                                "FN", credit);
                    } else { // PP
                        truth_vars->errtypes[truth_var_idx] = ERRTYPE_PP;
                        truth_vars->credit[truth_var_idx] = credit;
                        truth_vars->callq[truth_var_idx] = callq;
                        if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                truth_vars->refs[truth_var_idx].data(),
                                truth_vars->alts[truth_var_idx].data(), 
                                "PP", credit);
                    }
                }

                if (print) printf("SYNC @%s(%d,%d)=REF(%d,%d), ED %d->%d, QUERY %d-%d, TRUTH %d-%d\n",
                        (hi == ri) ? "REF" : "QUERY", qri, ti,
                        (hi == ri) ? qri : query_ref_ptrs[i][PTRS][qri],
                        truth_ref_ptrs[i][PTRS][ti], old_ed, new_ed, 
                        prev_query_var_ptr, query_var_ptr,
                        prev_truth_var_ptr, truth_var_ptr
                );

                prev_query_var_ptr = query_var_ptr;
                prev_truth_var_ptr = truth_var_ptr;
                new_ed = 0;
            }

            // update pointers and edit distance
            if (print) printf("%s %s (%d,%d) %s\n", 
                    sync[i][pidx] ? "*" : " ",
                    (hi == ri) ? "REF  " : "QUERY", qri, ti,
                    edits[i][pidx] ? "X" : "=");

            // update path pointer
            pidx--;
            if (pidx < 0) break;

            // update location
            qri = path[i][pidx].qri;
            ti = path[i][pidx].ti;
            hi = path[i][pidx].hi;
            new_ed += edits[i][pidx];

            // update next variant position
            query_var_pos = (query_var_ptr < query_beg_idx) ? -1 :
                query_vars->poss[query_var_ptr] - beg;
            truth_var_pos = (truth_var_ptr < truth_beg_idx) ? -1 :
                truth_vars->poss[truth_var_ptr] - beg;
        }
    }
}


/******************************************************************************/


void calc_edit_dist_aln(
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & offs,
        std::vector< std::vector< std::vector<int> > > & ptrs
        ) {

    // ALIGNMENT
    std::vector<std::string> query {query1, query1, query2, query2};
    std::vector<std::string> truth {truth1, truth2, truth1, truth2};
    std::vector<int> query_lens = 
            {int(query1.size()), int(query1.size()), int(query2.size()), int(query2.size())};
    std::vector<int> truth_lens = 
            {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};

    // for each combination of query and truth
    for(int i = 0; i < 4; i++) {

        int mat_len = query_lens[i] + truth_lens[i] - 1;
        offs[i].push_back(std::vector<int>(mat_len,-2));
        offs[i][0][query_lens[i]-1] = -1;
        ptrs[i].push_back(std::vector<int>(mat_len,PTR_NONE));
        bool done = false;
        while (true) {

            // EXTEND WAVEFRONT
            for (int d = 0; d < mat_len; d++) {
                int off = offs[i][s[i]][d];
                int diag = d + 1 - query_lens[i];

                // don't allow starting from untouched cells
                if (off == -2) continue;

                // check that it's within matrix
                if (diag + off + 1 < 0) continue;
                if (off > query_lens[i] - 1) continue;
                if (diag + off > truth_lens[i] - 1) continue;

                // extend
                while (off < query_lens[i] - 1 && 
                       diag + off < truth_lens[i] - 1) {
                    if (query[i][off+1] == truth[i][diag+off+1]) off++;
                    else break;
                }
                offs[i][s[i]][d] = off;

                // finish if done
                if (off == query_lens[i] - 1 && 
                    off + diag == truth_lens[i] - 1)
                { done = true; break; }

            }
            if (done) break;


            // NEXT WAVEFRONT
            // add wavefront, fill edge cells
            offs[i].push_back(std::vector<int>(mat_len, -2));
            // bottom left cells
            if (s[i]+1 == query_lens[i]-1)
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
}


/******************************************************************************/


editData alignment_wrapper(std::shared_ptr<superclusterData> clusterdata_ptr) {
    INFO(" ");
    INFO("Calculating precision/recall and distance metrics");

    // +2 since it's inclusive, but then also needs to include one quality higher
    // which doesn't contain any variants (to get draft reference edit dist)
    std::vector<int> all_qual_dists(g.max_qual+2, 0);
    editData edits;
    for (std::string ctg : clusterdata_ptr->contigs) {
        std::vector<int> ctg_qual_dists(g.max_qual+2,0);
        if (clusterdata_ptr->ctg_superclusters[ctg]->n)
            INFO("  Contig '%s'", ctg.data())

        // set superclusters pointer
        std::shared_ptr<ctgSuperclusters> sc = clusterdata_ptr->ctg_superclusters[ctg];

        // iterate over superclusters
        for(int sc_idx = 0; sc_idx < clusterdata_ptr->ctg_superclusters[ctg]->n; sc_idx++) {

            /////////////////////////////////////////////////////////////////////
            // PRECISION-RECALL: allow skipping called variants                  
            /////////////////////////////////////////////////////////////////////
            
            // set pointers between each hap (query1/2, truth1/2) and reference
            std::string query1 = "", ref_q1 = ""; 
            std::vector< std::vector<int> > query1_ref_ptrs, ref_query1_ptrs;
            generate_ptrs_strs(
                    query1, ref_q1, query1_ref_ptrs, ref_query1_ptrs, 
                    sc->ctg_variants[QUERY][HAP1], 
                    sc->superclusters[QUERY][HAP1][sc_idx],
                    sc->superclusters[QUERY][HAP1][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string query2 = "", ref_q2 = ""; 
            std::vector< std::vector<int> > query2_ref_ptrs, ref_query2_ptrs;
            generate_ptrs_strs(
                    query2, ref_q2, query2_ref_ptrs, ref_query2_ptrs, 
                    sc->ctg_variants[QUERY][HAP2],
                    sc->superclusters[QUERY][HAP2][sc_idx],
                    sc->superclusters[QUERY][HAP2][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth1 = "", ref_t1 = ""; 
            std::vector< std::vector<int> > truth1_ref_ptrs, ref_truth1_ptrs;
            generate_ptrs_strs(
                    truth1, ref_t1, truth1_ref_ptrs, ref_truth1_ptrs, 
                    sc->ctg_variants[TRUTH][HAP1],
                    sc->superclusters[TRUTH][HAP1][sc_idx],
                    sc->superclusters[TRUTH][HAP1][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth2 = "", ref_t2 = ""; 
            std::vector< std::vector<int> > truth2_ref_ptrs, ref_truth2_ptrs;
            generate_ptrs_strs(
                    truth2, ref_t2, truth2_ref_ptrs, ref_truth2_ptrs, 
                    sc->ctg_variants[TRUTH][HAP2],
                    sc->superclusters[TRUTH][HAP2][sc_idx],
                    sc->superclusters[TRUTH][HAP2][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );

            // calculate four forward-pass alignment edit dists
            // query1-truth2, query1-truth1, query2-truth1, query2-truth2
            std::vector<int> aln_score(4), aln_query_ref_end(4);
            std::vector< std::vector< std::vector<int> > > aln_ptrs;
            calc_prec_recall_aln(
                    query1, query2, truth1, truth2, ref_q1,
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs,
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    aln_score, aln_ptrs, aln_query_ref_end
            );

            // store optimal phasing for each supercluster
            // ORIG: query1-truth1 and query2-truth2
            // SWAP: query1-truth2 and query2-truth1
            int phase = store_phase(clusterdata_ptr, ctg, aln_score);

            // calculate paths from alignment
            std::vector< std::vector<idx1> > path;
            std::vector< std::vector<bool> > sync, edit;
            std::vector< std::vector< std::vector<int> > > path_ptrs, path_scores;
            calc_prec_recall_path(
                    clusterdata_ptr, sc_idx, ctg,
                    path, sync, edit,
                    aln_ptrs, path_ptrs, path_scores,
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs, 
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    aln_query_ref_end);

            // calculate precision/recall from paths
            calc_prec_recall(
                    clusterdata_ptr, sc_idx, ctg, ref_q1,
                    query1, query2, truth1, truth2,
                    path, sync, edit,
                    aln_ptrs, path_ptrs, 
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs,
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    aln_query_ref_end, phase, g.verbosity >= 2
            );


            /////////////////////////////////////////////////////////////////////
            // SMITH-WATERMAN DISTANCE: don't allow skipping called variants     
            /////////////////////////////////////////////////////////////////////
            
            // keep or swap truth haps based on previously decided phasing
            std::vector<std::string> truth(2);
            if (phase == PHASE_SWAP) {
                truth[HAP1] = truth2; truth[HAP2] = truth1;
            } else {
                truth[HAP1] = truth1; truth[HAP2] = truth2;
            }

            // phasing is known, add scores for each hap
            std::vector<int> qual_dists(g.max_qual+2, 0);
            for (int hap = 0; hap < HAPS; hap++) {

                // supercluster start/end indices
                int beg_idx = sc->ctg_variants[QUERY][hap]->clusters.size() ?
                    sc->ctg_variants[QUERY][hap]->clusters[
                        sc->superclusters[QUERY][hap][sc_idx]] : 0;
                int end_idx = sc->ctg_variants[QUERY][hap]->clusters.size() ?
                    sc->ctg_variants[QUERY][hap]->clusters[
                        sc->superclusters[QUERY][hap][sc_idx+1]] : 0;

                // calculate quality thresholds 
                // (where string would change when including/excluding variants)
                std::set<int> quals = {};
                for (int var_idx = beg_idx; var_idx < end_idx; var_idx++) {
                    quals.insert(sc->ctg_variants[QUERY][hap]->var_quals[var_idx]+1);
                }
                quals.insert(g.max_qual+2);

                // sweep through quality thresholds
                int prev_qual = 0;
                for (int qual : quals) {

                    // generate query string (only applying variants with Q>=qual)
                    std::string query = generate_str(
                            clusterdata_ptr->ref, 
                            sc->ctg_variants[QUERY][hap], 
                            ctg, beg_idx, end_idx,
                            sc->begs[sc_idx], sc->ends[sc_idx], 
                            prev_qual);

                    // align strings, backtrack, calculate distance
                    auto ptrs = sw_align(query, truth[hap], 
                            g.eval_sub, g.eval_open, g.eval_extend);
                    std::vector<int> cigar = sw_backtrack(query, truth[hap], ptrs);
                    int dist = count_dist(cigar);

                    // add distance for range of corresponding quals
                    for (int q = prev_qual; q < qual; q++) {
                        all_qual_dists[q] += dist;
                        ctg_qual_dists[q] += dist;
                        qual_dists[q] += dist;
                        edits.add_edits(ctg, sc->begs[sc_idx], hap, cigar, sc_idx, q);
                    }
                    prev_qual = qual;
                }
            }


            /////////////////////////////////////////////////////////////////////
            // DEBUG PRINTING                                                    
            /////////////////////////////////////////////////////////////////////
            if (g.verbosity >= 1) {
                // print cluster info
                printf("\n\nSupercluster: %d\n", sc_idx);
                for (int i = 0; i < CALLSETS*HAPS; i++) {
                    int callset = i >> 1;
                    int hap = i % 2;
                    int cluster_beg = sc->superclusters[callset][hap][sc_idx];
                    int cluster_end = sc->superclusters[callset][hap][sc_idx+1];
                    printf("%s%d: %d clusters (%d-%d)\n", 
                        callset_strs[callset].data(), hap+1,
                        cluster_end-cluster_beg,
                        cluster_beg, cluster_end);

                    for (int j = cluster_beg; j < cluster_end; j++) {
                        auto vars = sc->ctg_variants[callset][hap];
                        int variant_beg = vars->clusters[j];
                        int variant_end = vars->clusters[j+1];
                        printf("\tCluster %d: %d variants (%d-%d)\n", j, 
                            variant_end-variant_beg, variant_beg, variant_end);
                        for (int k = variant_beg; k < variant_end; k++) {
                            printf("\t\t%s %d\t%s\t%s\tQ=%f\n", ctg.data(), vars->poss[k], 
                            vars->refs[k].size() ?  vars->refs[k].data() : "_", 
                            vars->alts[k].size() ?  vars->alts[k].data() : "_",
                            vars->var_quals[k]);
                        }
                    }
                }

                printf("ORIG:      %s\n", ref_q1.data());
                printf("QUERY1:    %s\n", query1.data());
                printf("QUERY2:    %s\n", query2.data());
                printf("TRUTH1:    %s\n", truth1.data());
                printf("TRUTH2:    %s\n", truth2.data());
                printf("Edit Distance: %d\n", 
                    *std::min_element(qual_dists.begin(), qual_dists.end()));

            } // debug print

        } // each cluster
        if (clusterdata_ptr->ctg_superclusters[ctg]->n)
            INFO("    %d edits", 
                *std::min_element(ctg_qual_dists.begin(), ctg_qual_dists.end()));
    } // each contig
    INFO(" ");
    INFO("  Total edit distance: %d", 
                *std::min_element(all_qual_dists.begin(), all_qual_dists.end()));
    return edits;
}


/******************************************************************************/


/* Calculate the combined Smith-Waterman score for all variants within 
 * one or several adjacent clusters.
 */
int calc_vcf_sw_score(std::shared_ptr<ctgVariants> vars, 
        int clust_beg_idx, int clust_end_idx, int sub, int open, int extend) {
    int score = 0;
    for (int var_idx = vars->clusters[clust_beg_idx]; 
            var_idx < vars->clusters[clust_end_idx]; var_idx++) {
        switch (vars->types[var_idx]) {
            case TYPE_SUB:
                score += sub;
                break;
            case TYPE_INS:
                score += open;
                score += extend * vars->alts[var_idx].size();
                break;
            case TYPE_DEL:
                score += open;
                score += extend * vars->refs[var_idx].size();
                break;
            default:
                ERROR("Unexpected variant type in calc_vcf_sw_score()");
        }
    }
    return score;
}


int calc_cig_sw_score(const std::vector<int> & cigar, 
        int sub, int open, int extend) {
    int score = 0;
    int prev = PTR_DIAG;
    for (int i = 0; i < int(cigar.size()); i++) {
        switch (cigar[i]) {
            case PTR_DIAG:
                i++; // increment by 2
                prev = PTR_DIAG;
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
                ERROR("Unexpected variant type in calc_cig_sw_score()");
        }
    }
    return score;
}


/******************************************************************************/


/* Perform Djikstra Smith-Waterman alignment of two strings, returning 
 * the farthest-reaching reference index of lesser or equal score to that provided.
 * Strings are generated by applying variants to draft ref, and skipping 
 * variants is not allowed. Neither is a diagonal ref transition.
 */
int sw_max_reach(const std::string & query, const std::string & ref, 
        const std::vector< std::vector<int> > & query_ref_ptrs,
        const std::vector< std::vector<int> > & ref_query_ptrs,
        int sub, int open, int extend,
        int score, bool reverse /*= false*/, int ref_section /*= -1*/) {
    
    int ref_len = ref.size();
    int query_len = query.size();
    int sub_wave = 0, open_wave = 0, extend_wave = 0, end_wave = 0;
        
    // set first wavefront
    std::queue<idx2> queue;
    queue.push({MAT_SUB, 0, 0});
    std::unordered_set<idx2> done;
    done.insert({MAT_SUB, 0, 0});
    std::vector< std::unordered_set<idx2> > waves(score+1);

    // find furthest-reaching alignment with lesser or equal score
    for (int s = 0; s <= score; s++) {
        if (g.verbosity >= 3)
            printf("s = %d\n", s);

        // EXTEND WAVEFRONT (stay at same score)
        while (!queue.empty()) {
            idx2 x = queue.front(); queue.pop();
            if (g.verbosity >= 3)
                printf("  x = (%c, %d, %d)\n", 
                    std::string("SID")[x.mi], x.qi, x.ri);
            waves[s].insert(x);

            // allow non-main-diagonal match
            bool main_diag = x.ri >= ref_section &&
                    !(query_ref_ptrs[FLAGS][x.qi] & PTR_VARIANT) && // neither prev/curr are variants
                    !(ref_query_ptrs[FLAGS][x.ri] & PTR_VARIANT) &&
                    !(query_ref_ptrs[FLAGS][x.qi+1] & PTR_VARIANT) &&
                    !(ref_query_ptrs[FLAGS][x.ri+1] & PTR_VARIANT) &&
                    query_ref_ptrs[PTRS][x.qi] == x.ri && // point to one another (main diag)
                    ref_query_ptrs[PTRS][x.ri] == x.qi &&
                    query_ref_ptrs[PTRS][x.qi+1] == x.ri+1 &&
                    ref_query_ptrs[PTRS][x.ri+1] == x.qi+1;
            if (x.qi+1 < query_len && x.ri+1 < ref_len && 
                    !main_diag &&
                    query[x.qi+1] == ref[x.ri+1] && 
                    !contains(done, {x.mi, x.qi+1, x.ri+1})) {
                idx2 next(x.mi, x.qi+1, x.ri+1);
                queue.push(next);
                done.insert(next);
            }

            // allow exiting D/I state (no penalty for forward only)
            if (!reverse && (x.mi == MAT_INS || x.mi == MAT_DEL) && 
                    !contains(done, {MAT_SUB, x.qi, x.ri})) {
                queue.push({MAT_SUB, x.qi, x.ri});
                done.insert({MAT_SUB, x.qi, x.ri});
            }
        }

        // exit if we're done aligning
        if (s == score || contains(done, {MAT_SUB, query_len-1, ref_len-1})) 
            break;

        // NEXT WAVEFRONT (increase score by one)
        
        // SUB transition (only from MAT_SUB)
        sub_wave = s + 1 - sub;
        if (sub_wave >= 0 && sub_wave <= score) {
            for (idx2 x : waves[sub_wave]) {
                if (x.mi == MAT_SUB && x.qi+1 < query_len && x.ri+1 < ref_len) {
                    idx2 next({x.mi, x.qi+1, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }
            }
        }

        // INS/DEL opening (only from MAT_SUB)
        if (reverse) // g.open penalty when leaving INDEL
            open_wave = s + 1 - extend;
        else // g.open penalty when entering INDEL
            open_wave = s + 1 - (open+extend);
        if (open_wave >= 0 && open_wave <= score) {
            for(idx2 x : waves[open_wave]) {

                // INS opening
                if (x.mi == MAT_SUB && x.qi+1 < query_len) {
                    idx2 next({MAT_INS, x.qi+1, x.ri});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }

                // DEL opening
                if (x.mi == MAT_SUB && x.ri+1 < ref_len) {
                    idx2 next({MAT_DEL, x.qi, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }
            }
        }

        // reverse alignment, penalize leaving INDEL
        end_wave = s + 1 - open;
        if (reverse && end_wave >= 0 && end_wave <= score) {
            for (idx2 x : waves[end_wave]) {
                if ((x.mi == MAT_INS || x.mi == MAT_DEL) && 
                        !contains(done, {MAT_SUB, x.qi, x.ri})) {
                    queue.push({MAT_SUB, x.qi, x.ri});
                    done.insert({MAT_SUB, x.qi, x.ri});
                }
            }
        }

        // INS/DEL extension
        extend_wave = s + 1 - extend;
        if (extend_wave >= 0 && extend_wave <= score) {
            for(idx2 x : waves[extend_wave]) {

                // INS extending (only from MAT_INS)
                if (x.mi == MAT_INS && x.qi+1 < query_len) {
                    idx2 next({x.mi, x.qi+1, x.ri});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }

                // DEL extending (only from MAT_DEL)
                if (x.mi == MAT_DEL && x.ri+1 < ref_len) {
                    idx2 next({x.mi, x.qi, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }
            }
        }
    } // end reach

    // farthest reach isn't necessarily at max score
    int min_wave = std::max(0, std::min(std::min(
                    open_wave, extend_wave), sub_wave));
    // get max reach
    int max_reach = 0;
    for (int w = min_wave; w <= score; w++) {
        for (idx2 x : waves[w]) {
            max_reach = std::max(max_reach, x.ri);
        }
    }
    return max_reach;
}


/******************************************************************************/


std::unique_ptr<variantData> sw_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref_fasta, 
        int sub, int open, int extend) {
    INFO(" ");
    INFO("Realigning VCF '%s'", vcf->filename.data());

    // copy vcf header data over to results vcf
    std::unique_ptr<variantData> results(new variantData());
    results->set_header(vcf);

    // iterate over each contig haplotype
    for (int hap = 0; hap < 2; hap++) {
        for (auto itr = vcf->ctg_variants[hap].begin(); 
                itr != vcf->ctg_variants[hap].end(); itr++) {
            std::string ctg = itr->first;
            std::shared_ptr<ctgVariants> vars = itr->second;
            if (vars->poss.size() == 0) continue;
            INFO("  Haplotype %d Contig %s", hap+1, ctg.data());

            // realign each cluster of variants
            for (int cluster = 0; cluster < int(vars->clusters.size()-1); cluster++) {
                int beg_idx = vars->clusters[cluster];
                int end_idx = vars->clusters[cluster+1];
                int beg = vars->poss[beg_idx]-1;
                int end = vars->poss[end_idx-1] + vars->rlens[end_idx-1]+1;

                // variant qual is minimum in cluster
                float qual = g.max_qual;
                for (int i = beg_idx; i < end_idx; i++) {
                    qual = std::min(qual, vars->var_quals[i]);
                }

                // generate strings
                std::string query = 
                    generate_str(ref_fasta, vars, ctg, beg_idx, end_idx, beg, end);
                std::string ref = ref_fasta->fasta.at(ctg).substr(beg, end-beg);
                
                // perform alignment
                if (g.verbosity >= 2) {
                    printf("REF:   %s\n", ref.data());
                    printf("QUERY: %s\n", query.data());
                }
                std::unordered_map<idx2, idx2> ptrs = sw_align(query, ref,
                        sub, open, extend);
                
                // backtrack
                std::vector<int> cigar = sw_backtrack(query, ref, ptrs);
                if (g.verbosity >= 2)
                    print_cigar(cigar);

                // compare distances
                if (g.verbosity >= 1) {
                    int new_score = calc_cig_sw_score(cigar, sub, open, extend);
                    int old_score = calc_vcf_sw_score(vars, cluster, cluster+1, 
                            sub, open, extend);
                    if (new_score < old_score) {
                        printf("\n\tCluster %d: %d variants (%d-%d)\n", 
                            cluster, end_idx-beg_idx, beg_idx, end_idx);
                        for (int i = beg_idx; i < end_idx; i++) {
                            printf("\t\t%s %d\t%s\t%s\tQ=%f\n", 
                                ctg.data(), vars->poss[i], 
                                vars->refs[i].size() ?  vars->refs[i].data() : "_", 
                                vars->alts[i].size() ?  vars->alts[i].data() : "_",
                                vars->var_quals[i]);
                        }
                        printf("Old score: %d\n", old_score);
                        printf("New score: %d\n", new_score);
                    }
                }
                
                // save resulting variants
                results->add_variants(cigar, hap, beg, ctg, query, ref, qual);

            } // cluster
        } // contig
    } // hap

    return results;
}


/******************************************************************************/


/* Perform global Djikstra Smith-Waterman alignment of two strings. 
 * Results are stored in `ptrs`.
 * */
std::unordered_map<idx2, idx2> sw_align(
        const std::string & query, 
        const std::string & ref, 
        int sub, int open, int extend) { 
    
    std::unordered_map<idx2, idx2> ptrs;
    ptrs[idx2(0,0,0)] = idx2(0,-1,-1);
    int ref_len = ref.size();
    int query_len = query.size();
    int sub_wave = 0, open_wave = 0, extend_wave = 0;
        
    // set first wavefront
    std::queue<idx2> queue;
    queue.push({MAT_SUB, 0, 0});
    std::unordered_set<idx2> done;
    done.insert({MAT_SUB, 0, 0});
    std::vector< std::set<idx2> > waves;

    // perform global alignment
    int s = 0;
    while(true) {
        waves.push_back(std::set<idx2>());
        if (g.verbosity >= 3)
            printf("s = %d\n", s);

        // EXTEND WAVEFRONT (stay at same score)
        while (!queue.empty()) {
            idx2 x = queue.front(); queue.pop();
            idx2 prev = ptrs.find(x)->second;
            if (g.verbosity >= 3)
                printf("  x = (%c, %d, %d) -> (%c, %d, %d)\n", 
                    std::string("SID")[x.mi], x.qi, x.ri,
                    std::string("SID")[prev.mi], prev.qi, prev.ri);
            waves[s].insert(x);

            // allow match
            if (x.mi == MAT_SUB &&
                    x.qi+1 < query_len && x.ri+1 < ref_len &&
                    query[x.qi+1] == ref[x.ri+1] && 
                    !contains(done, {x.mi, x.qi+1, x.ri+1})) {
                idx2 next(x.mi, x.qi+1, x.ri+1);
                queue.push(next); done.insert(next); ptrs[next] = x;
            }

            // allow exiting D/I state freely
            if ((x.mi == MAT_INS || x.mi == MAT_DEL) && 
                    !contains(done, {MAT_SUB, x.qi, x.ri})) {
                idx2 next(MAT_SUB, x.qi, x.ri);
                queue.push(next); done.insert(next); ptrs[next] = x;
            }
        }

        // exit if we're done aligning
        if (contains(done, {MAT_SUB, query_len-1, ref_len-1})) 
            break;

        // NEXT WAVEFRONT (increase score by one)
        
        // SUB transition (only from MAT_SUB)
        sub_wave = s + 1 - sub;
        if (sub_wave >= 0) {
            for (idx2 x : waves[sub_wave]) {
                if (x.mi == MAT_SUB && x.qi+1 < query_len && x.ri+1 < ref_len) {
                    idx2 next(x.mi, x.qi+1, x.ri+1);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;

                        // continue extending while match
                        idx2 prev(next);
                        next.qi++; next.ri++;
                        while (next.qi < query_len && next.ri < ref_len &&
                                query[next.qi] == ref[next.ri] &&
                                !contains(done, next)) {
                            queue.push(next); done.insert(next); ptrs[next] = prev;

                            // update
                            prev = next;
                            next.qi++; next.ri++;
                        }
                    }
                }
            }
        }

        // INS/DEL opening (only from MAT_SUB)
        open_wave = s + 1 - (open+extend);
        if (open_wave >= 0) {
            for(idx2 x : waves[open_wave]) {

                // INS opening
                if (x.mi == MAT_SUB && x.qi+1 < query_len) {
                    idx2 next(MAT_INS, x.qi+1, x.ri);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;

                        // shift back to SUB_MAT immediately
                        idx2 prev(next);
                        next.mi = MAT_SUB;
                        if (!contains(done, next)) {
                            queue.push(next); done.insert(next); ptrs[next] = prev;

                            // continue extending while match
                            prev = next;
                            next.qi++; next.ri++;
                            while (next.qi < query_len && next.ri < ref_len &&
                                    query[next.qi] == ref[next.ri] &&
                                    !contains(done, next)) {
                                queue.push(next); done.insert(next); ptrs[next] = prev;

                                // update
                                prev = next;
                                next.qi++; next.ri++;
                            }
                        }
                    }
                }

                // DEL opening
                if (x.mi == MAT_SUB && x.ri+1 < ref_len) {
                    // add to queue if in range and unvisited
                    idx2 next(MAT_DEL, x.qi, x.ri+1);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;

                        // shift back to SUB_MAT immediately
                        idx2 prev(next);
                        next.mi = MAT_SUB;
                        if (!contains(done, next)) {
                            queue.push(next); done.insert(next); ptrs[next] = prev;

                            // continue extending while match
                            prev = next;
                            next.qi++; next.ri++;
                            while (next.qi < query_len && next.ri < ref_len &&
                                    query[next.qi] == ref[next.ri] &&
                                    !contains(done, next)) {
                                queue.push(next); done.insert(next); ptrs[next] = prev;

                                // update
                                prev = next;
                                next.qi++; next.ri++;
                            }
                        }
                    }
                }
            }
        }

        // INS/DEL extension
        extend_wave = s + 1 - extend;
        if (extend_wave >= 0) {
            for(idx2 x : waves[extend_wave]) {

                // INS extending (only from MAT_INS)
                if (x.mi == MAT_INS && x.qi+1 < query_len) {
                    // add to queue if in range and unvisited
                    idx2 next(MAT_INS, x.qi+1, x.ri);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;

                        // shift back to SUB_MAT immediately
                        idx2 prev(next);
                        next.mi = MAT_SUB;
                        if (!contains(done, next)) {
                            queue.push(next); done.insert(next); ptrs[next] = prev;

                            // continue extending while match
                            prev = next;
                            next.qi++; next.ri++;
                            while (next.qi < query_len && next.ri < ref_len &&
                                    query[next.qi] == ref[next.ri] &&
                                    !contains(done, next)) {
                                queue.push(next); done.insert(next); ptrs[next] = prev;

                                // update
                                prev = next;
                                next.qi++; next.ri++;
                            }
                        }
                    }
                }

                // DEL extending (only from MAT_DEL)
                if (x.mi == MAT_DEL && x.ri+1 < ref_len) {
                    // add to queue if in range and unvisited
                    idx2 next(MAT_DEL, x.qi, x.ri+1);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;

                        // shift back to SUB_MAT immediately
                        idx2 prev(next);
                        next.mi = MAT_SUB;
                        if (!contains(done, next)) {
                            queue.push(next); done.insert(next); ptrs[next] = prev;

                            // continue extending while match
                            prev = next;
                            next.qi++; next.ri++;
                            while (next.qi < query_len && next.ri < ref_len &&
                                    query[next.qi] == ref[next.ri] &&
                                    !contains(done, next)) {
                                queue.push(next); done.insert(next); ptrs[next] = prev;

                                // update
                                prev = next;
                                next.qi++; next.ri++;
                            }
                        }
                    }
                }
            }
        }

        s++;
    } // end alignment
    return ptrs;
}


/******************************************************************************/


int count_dist(const std::vector<int> & cigar) {
    int cigar_ptr = 0;
    int dist = 0;
    while (cigar_ptr < int(cigar.size())) {
        switch (cigar[cigar_ptr]) {
        case PTR_DIAG:
            cigar_ptr += 2;
            break;
        case PTR_SUB:
            dist++;
            cigar_ptr += 2;
            break;
        case PTR_INS:
        case PTR_DEL:
            dist++;
            cigar_ptr++;
            break;
        }
    }
    return dist;
}


/******************************************************************************/


/* After alignment, backtrack using pointers and save cigar. */
std::vector<int> sw_backtrack(
        const std::string & query,
        const std::string & ref,
        const std::unordered_map<idx2, idx2> & ptrs) {

    // initialize backtrack position and cigar string
    std::vector<int> cigar(query.size() + ref.size(), 0);
    idx2 pos(MAT_SUB, query.size()-1, ref.size()-1);
    int cigar_ptr = cigar.size()-1;

    if (g.verbosity >= 2) {
        printf("Backtrack:\n");
    }
    while (pos.qi >= 0 || pos.ri >= 0) {
        if (g.verbosity >= 2) {
            printf("(%c, %d, %d)\n",
                    std::string("SID")[pos.mi], pos.qi, pos.ri);
        }
        auto next_itr = ptrs.find(pos);
        if (next_itr == ptrs.end())
            ERROR("(%c, %d, %d) not found in ptrs during backtrack.",
                    std::string("SID")[pos.mi], pos.qi, pos.ri);
        idx2 next = ptrs.find(pos)->second;
        if (next.mi == pos.mi) { // stay in same matrix
            switch (pos.mi) {
                case MAT_SUB:
                    if (query[pos.qi] == ref[pos.ri]) { // match
                        cigar[cigar_ptr--] = PTR_DIAG;
                        cigar[cigar_ptr--] = PTR_DIAG;
                    } else { // SUB
                        cigar[cigar_ptr--] = PTR_SUB;
                        cigar[cigar_ptr--] = PTR_SUB;
                    }
                    break;
                case MAT_INS: // extend INS
                    cigar[cigar_ptr--] = PTR_INS;
                    break;
                case MAT_DEL: // extend DEL
                    cigar[cigar_ptr--] = PTR_DEL;
                    break;
                default:
                    ERROR("Unexpected MAT type in sw_backtrack()");
            }
        } else { // change matrices
            switch (pos.mi) {
                case MAT_SUB: // do nothing, left INDEL
                    break;
                case MAT_INS: // enter INS
                    cigar[cigar_ptr--] = PTR_INS;
                    break;
                case MAT_DEL: // enter DEL
                    cigar[cigar_ptr--] = PTR_DEL;
                    break;
                default:
                    ERROR("Unexpected MAT type in sw_backtrack()");
            }
        }
        pos = next;
    }
    return cigar;
}
