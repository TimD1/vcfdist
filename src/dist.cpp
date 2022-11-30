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
                    // TODO: VALGRIND: conditional jump or move depends on uninitialized value
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
                        int alt_idx = cell.ci;
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
                        int alt_idx = cell.ci;
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
void reverse_ptrs_strs(std::string & query, std::string & ref,
        std::vector<int> & query_ptrs, std::vector<int> & ref_ptrs) {

    // reverse strings and vectors
    std::reverse(query.begin(), query.end());
    std::reverse(ref.begin(), ref.end());
    std::reverse(query_ptrs.begin(), query_ptrs.end());
    std::reverse(ref_ptrs.begin(), ref_ptrs.end());

    // update pointers
    for(int & ptr : query_ptrs) {
        if (ptr >= 0) { ptr = ref.size()-1 - ptr; }
    }
    for(int & ptr : ref_ptrs) {
        if (ptr >= 0) { ptr = query.size()-1 - ptr; }
    }
    return;
}


/******************************************************************************/


/* Generate the new sequence by applying variants to the reference. */
std::string generate_str(
        std::shared_ptr<fastaData> ref, 
        std::shared_ptr<ctgVariants> vars, std::string ctg,
        int var_idx, int end_idx, int beg_pos, int end_pos, int min_qual=0) {

    std::string str = "";
    for (int ref_pos = beg_pos; ref_pos < end_pos; ) {

        // VARIANT
        if (ref_pos == vars->poss[var_idx]  && var_idx < end_idx) {
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
                    case TYPE_GRP:
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
                ERROR("Contig %s not present in reference FASTA",
                        ctg.data());
                exit(1);
            }
        }
    }
    return str;
}


/******************************************************************************/


/* For each position on hap1 and hap2, it generates pointers from one hap to the
 * other if both match the reference, or else -1. This method also generates the
 * strings for alignment, and the pretty-printing strings for debug output.
 */
void generate_ptrs_strs(
        std::string & hap1, std::string & hap2,          // actual strings 
        std::vector<int> & hap1_ptrs, std::vector<int> & hap2_ptrs,
        std::shared_ptr<ctgVariants> hap1_vars, std::shared_ptr<ctgVariants> hap2_vars,
        size_t hap1_clust_beg_idx, size_t hap2_clust_beg_idx,
        size_t hap1_clust_end_idx, size_t hap2_clust_end_idx,
        int beg_pos, int end_pos, std::shared_ptr<fastaData> ref, std::string ctg
        ) {

    // generate hap1 and hap2 strings and pointers
    int hap1_var_idx = hap1_vars->clusters[hap1_clust_beg_idx];
    int hap2_var_idx = hap2_vars->clusters[hap2_clust_beg_idx];
    int hap1_end_idx = hap1_vars->clusters[hap1_clust_end_idx];
    int hap2_end_idx = hap2_vars->clusters[hap2_clust_end_idx];
    for (int hap1_ref_pos = beg_pos, hap2_ref_pos = beg_pos; 
            hap1_ref_pos < end_pos || hap2_ref_pos < end_pos; ) {

        // CONSIDER TRUTH1 ONLY, PRIOR REFERENCE POSITION
        if (hap1_ref_pos < hap2_ref_pos) {
            if (hap1_ref_pos == hap1_vars->poss[hap1_var_idx] && 
                    hap1_var_idx < hap1_end_idx) { // in hap1 variant
                switch (hap1_vars->types[hap1_var_idx]) {
                    case TYPE_INS:
                        hap1 += hap1_vars->alts[hap1_var_idx];
                        hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                        break;
                    case TYPE_DEL:
                        hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                        break;
                    case TYPE_SUB:
                        hap1 += hap1_vars->alts[hap1_var_idx];
                        hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                        hap1_ref_pos++;
                        break;
                    case TYPE_GRP:
                        hap1 += hap1_vars->alts[hap1_var_idx];
                        hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                        hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                        break;
                }
                hap1_var_idx++; // next variant
            } else { // no hap1 variant, in hap2 variant
                try {
                    hap1 += ref->fasta.at(ctg)[hap1_ref_pos];
                    hap1_ptrs.push_back(-1);
                    hap1_ref_pos++;
                } catch (const std::out_of_range & e) {
                    ERROR("Contig %s not present in reference FASTA",
                            ctg.data());
                    exit(1);
                }
            }
        }

        // CONSIDER TRUTH2 ONLY, PRIOR REFERENCE POSITION
        else if (hap2_ref_pos < hap1_ref_pos) {
            if (hap2_ref_pos == hap2_vars->poss[hap2_var_idx] &&
                    hap2_var_idx < hap2_end_idx) { // in hap2 variant
                switch (hap2_vars->types[hap2_var_idx]) {
                    case TYPE_INS:
                        hap2 += hap2_vars->alts[hap2_var_idx];
                        hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                        break;
                    case TYPE_DEL:
                        hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                        break;
                    case TYPE_SUB:
                        hap2 += hap2_vars->alts[hap2_var_idx];
                        hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                        hap2_ref_pos++;
                        break;
                    case TYPE_GRP:
                        hap2 += hap2_vars->alts[hap2_var_idx];
                        hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                        hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                        break;
                }
                hap2_var_idx++; // next variant
            } else { // match
                try {
                    hap2 += ref->fasta.at(ctg)[hap2_ref_pos];
                    hap2_ptrs.push_back(-1);
                    hap2_ref_pos++;
                } catch (const std::out_of_range & e) {
                    ERROR("Contig %s not present in reference FASTA",
                            ctg.data());
                    exit(1);
                }
            }
        }

        // REFERENCE POSITIONS MATCH! POTENTIAL TRANSITIONS
        else {
            bool hap1_var = false;
            if (hap1_ref_pos == hap1_vars->poss[hap1_var_idx] && 
                    hap1_var_idx < hap1_end_idx) { // in hap1 variant
                hap1_var = true;
                switch (hap1_vars->types[hap1_var_idx]) {
                    case TYPE_INS:
                        hap1 += hap1_vars->alts[hap1_var_idx];
                        hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                        break;
                    case TYPE_DEL:
                        hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                        break;
                    case TYPE_SUB:
                        hap1 += hap1_vars->alts[hap1_var_idx];
                        hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                        hap1_ref_pos++;
                        break;
                    case TYPE_GRP:
                        hap1 += hap1_vars->alts[hap1_var_idx];
                        hap1_ptrs.insert(hap1_ptrs.end(), hap1_vars->alts[hap1_var_idx].size(), -1);
                        hap1_ref_pos += hap1_vars->refs[hap1_var_idx].size();
                        break;
                }
                hap1_var_idx++; // next variant
            } 

            bool hap2_var = false;
            if (hap2_ref_pos == hap2_vars->poss[hap2_var_idx] && 
                    hap2_var_idx < hap2_end_idx) { // in hap2 variant
                hap2_var = true;
                switch (hap2_vars->types[hap2_var_idx]) {
                    case TYPE_INS:
                        hap2 += hap2_vars->alts[hap2_var_idx];
                        hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                        break;
                    case TYPE_DEL:
                        hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                        break;
                    case TYPE_SUB:
                        hap2 += hap2_vars->alts[hap2_var_idx];
                        hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                        hap2_ref_pos++;
                        break;
                    case TYPE_GRP:
                        hap2 += hap2_vars->alts[hap2_var_idx];
                        hap2_ptrs.insert(hap2_ptrs.end(), hap2_vars->alts[hap2_var_idx].size(), -1);
                        hap2_ref_pos += hap2_vars->refs[hap2_var_idx].size();
                        break;
                }
                hap2_var_idx++; // next variant

            } 

            // BOTH MATCH REFERENCE, ADD POINTERS FOR TRANSITIONS
            if (!hap1_var && !hap2_var) { // add pointers
                try {

                    // find next position w/o ref match
                    int hap1_next_var = (hap1_var_idx < hap1_end_idx) ?
                        hap1_vars->poss[hap1_var_idx] : end_pos;
                    int hap2_next_var = (hap2_var_idx < hap2_end_idx) ?
                        hap2_vars->poss[hap2_var_idx] : end_pos;
                    int ref_end = std::min(
                            std::min(hap1_next_var, hap2_next_var), end_pos);

                    // add pointers
                    std::vector<int> new_hap1_ptrs(ref_end-hap1_ref_pos);
                    for(size_t i = 0; i < new_hap1_ptrs.size(); i++)
                        new_hap1_ptrs[i] = hap2.size() + i;
                    hap1_ptrs.insert(hap1_ptrs.end(), 
                            new_hap1_ptrs.begin(), new_hap1_ptrs.end());
                    std::vector<int> new_hap2_ptrs(ref_end-hap2_ref_pos);
                    for(size_t i = 0; i < new_hap2_ptrs.size(); i++)
                        new_hap2_ptrs[i] = hap1.size() + i;
                    hap2_ptrs.insert(hap2_ptrs.end(), 
                            new_hap2_ptrs.begin(), new_hap2_ptrs.end());

                    // add sequence, update positions
                    hap1 += ref->fasta.at(ctg).substr(hap1_ref_pos, ref_end-hap1_ref_pos);
                    hap2 += ref->fasta.at(ctg).substr(hap2_ref_pos, ref_end-hap2_ref_pos);
                    hap1_ref_pos = ref_end;
                    hap2_ref_pos = ref_end;

                } catch (const std::out_of_range & e) {
                    ERROR("Contig '%s' not present in reference FASTA",
                            ctg.data());
                    exit(1);
                }
            }
        }
    }
}


/******************************************************************************/


void calc_prec_recall_aln(
        std::string query1, std::string query2,
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & ptrs,
        std::vector<int> & pr_query_ref_end
        ) {
    
    // set loop variables
    int ref_len = ref.size();
    std::vector<std::string> query {query1, query1, query2, query2};
    std::vector<std::string> truth {truth1, truth2, truth1, truth2};
    std::vector< std::vector<int> > query_ref_ptrs {
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector<int> > ref_query_ptrs {
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector<int> query_lens = 
            {int(query1.size()), int(query1.size()), int(query2.size()), int(query2.size())};
    std::vector<int> truth_lens = 
            {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};

    std::vector< std::vector< std::vector<bool> > > done;

    // for each combination of query and truth
    for (int i = 0; i < 4; i++) {
        int ci = 2*i + QUERY; // query index (offs and ptrs)
        int ri = 2*i + REF;   // ref index   (offs and ptrs)

        // init full pointer/done matrices
        ptrs.push_back(std::vector< std::vector<int> >(query_lens[i], 
                std::vector<int>(truth_lens[i], PTR_NONE)));
        ptrs.push_back(std::vector< std::vector<int> >(ref_len, 
                std::vector<int>(truth_lens[i], PTR_NONE)));
        ptrs[ci][0][0] = PTR_DIAG;

        done.push_back(std::vector< std::vector<bool> >(query_lens[i], 
                std::vector<bool>(truth_lens[i], false)));
        done.push_back(std::vector< std::vector<bool> >(ref_len, 
                std::vector<bool>(truth_lens[i], false)));
        
        // set first wavefront
        std::queue<idx1> queue;
        queue.push({ci, 0, 0});
        ptrs[ci][0][0] |= PTR_DIAG;
        done[ci][0][0] = true;
        queue.push({ri, 0, 0});
        ptrs[ri][0][0] |= PTR_DIAG;
        done[ri][0][0] = true;

        // continue looping until full alignment found
        std::unordered_set<idx1> done_this_wave;
        std::unordered_set<idx1> this_wave;
        while (true) {
            /* printf("s = %d\n", s[i]); */
            if (queue.empty()) ERROR("Empty queue in 'prec_recall_aln()'.");

            // EXTEND WAVEFRONT (stay at same score)
            while (!queue.empty()) {
                idx1 x = queue.front(); queue.pop();
                /* printf("  x = (%d, %d, %d)\n", x.hi, x.cri, x.ti); */
                this_wave.insert(x);
                if (x.hi == ci) { // == QUERY
                    // allow match
                    if (x.cri+1 < query_lens[i] && x.ti+1 < truth_lens[i] &&
                            query[i][x.cri+1] == truth[i][x.ti+1]) {
                        if (!done[x.hi][x.cri+1][x.ti+1] && 
                                !contains(done_this_wave, idx1(x.hi, x.cri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.cri+1, x.ti+1));
                            done_this_wave.insert(idx1(x.hi, x.cri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.cri+1][x.ti+1]) {
                            ptrs[x.hi][x.cri+1][x.ti+1] |= PTR_DIAG;
                        }
                    }
                    // allow phase swap
                    if (query_ref_ptrs[i][x.cri] >= 0) {
                        if (!done[ri][query_ref_ptrs[i][x.cri]][x.ti] &&
                                !contains(done_this_wave, idx1(ri, query_ref_ptrs[i][x.cri], x.ti))) {
                            queue.push(idx1(ri, 
                                    query_ref_ptrs[i][x.cri], 
                                    x.ti));
                            done_this_wave.insert(idx1(ri, query_ref_ptrs[i][x.cri], x.ti));
                        }
                        if (!done[ri][query_ref_ptrs[i][x.cri]][x.ti]) {
                            ptrs[ri][query_ref_ptrs[i][x.cri]][x.ti] |= PTR_SWAP;
                        }
                    }
                } else { // x.hi == ri == REF
                    // allow match
                    if (x.cri+1 < ref_len && x.ti+1  < truth_lens[i] &&
                            ref[x.cri+1] == truth[i][x.ti+1]) {
                        if (!done[x.hi][x.cri+1][x.ti+1] &&
                                !contains(done_this_wave, idx1(x.hi, x.cri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.cri+1, x.ti+1));
                            done_this_wave.insert(idx1(x.hi, x.cri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.cri+1][x.ti+1]) {
                            ptrs[x.hi][x.cri+1][x.ti+1] |= PTR_DIAG;
                        }
                    }
                    // allow phase swap
                    if (ref_query_ptrs[i][x.cri] >= 0) {
                        if (!done[ci][ref_query_ptrs[i][x.cri]][x.ti] &&
                                !contains(this_wave, idx1(ci, ref_query_ptrs[i][x.cri], x.ti))) {
                            queue.push(idx1(ci, 
                                        ref_query_ptrs[i][x.cri], 
                                        x.ti));
                            done_this_wave.insert(idx1(ci, ref_query_ptrs[i][x.cri], x.ti));
                        }
                        if (!done[ci][ref_query_ptrs[i][x.cri]][x.ti]) {
                            ptrs[ci][ref_query_ptrs[i][x.cri]][x.ti] |= PTR_SWAP;
                        }
                    }
                }
            }

            // mark all cells visited this wave as done
            for (auto x : done_this_wave) { done[x.hi][x.cri][x.ti] = true; }
            done_this_wave.clear();

            // exit if we're done aligning
            if (done[ci][query_lens[i]-1][truth_lens[i]-1] ||
                done[ri][ref_len-1][truth_lens[i]-1]) break;


            // NEXT WAVEFRONT (increase score by one)
            for (auto x : this_wave) {
                if (x.hi == ci) { // QUERY
                    if (x.cri+1 < query_lens[i]) { // INS
                        if (!done[x.hi][x.cri+1][x.ti] && 
                                !contains(done_this_wave,
                                    idx1(x.hi, x.cri+1, x.ti))) {
                            queue.push(idx1(x.hi, x.cri+1, x.ti));
                            done_this_wave.insert(idx1(x.hi, x.cri+1, x.ti));
                        }
                        if (!done[x.hi][x.cri+1][x.ti])
                            ptrs[x.hi][x.cri+1][x.ti] |= PTR_UP;
                    }
                    if (x.ti+1 < truth_lens[i]) { // DEL
                        if (!done[x.hi][x.cri][x.ti+1] &&
                                !contains(done_this_wave,
                                    idx1(x.hi, x.cri, x.ti+1))) {
                            queue.push(idx1(x.hi, x.cri, x.ti+1));
                            done_this_wave.insert(idx1(x.hi, x.cri, x.ti+1));
                        }
                        if (!done[x.hi][x.cri][x.ti+1])
                            ptrs[x.hi][x.cri][x.ti+1] |= PTR_LEFT;
                    }
                    if (x.cri+1 < query_lens[i] && x.ti+1 < truth_lens[i]) { // SUB
                        if (!done[x.hi][x.cri+1][x.ti+1] &&
                                !contains(done_this_wave,
                                    idx1(x.hi, x.cri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.cri+1, x.ti+1));
                            done_this_wave.insert(idx1(x.hi, x.cri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.cri+1][x.ti+1])
                            ptrs[x.hi][x.cri+1][x.ti+1] |= PTR_SUB;
                    }
                } else { //x.hi == REF
                    if (x.cri+1 < ref_len) { // INS
                        if (!done[x.hi][x.cri+1][x.ti] &&
                                !contains(done_this_wave,
                                    idx1(x.hi, x.cri+1, x.ti))) {
                            queue.push(idx1(x.hi, x.cri+1, x.ti));
                            done_this_wave.insert(idx1(x.hi, x.cri+1, x.ti));
                        }
                        if (!done[x.hi][x.cri+1][x.ti])
                            ptrs[x.hi][x.cri+1][x.ti] |= PTR_UP;
                    }
                    if (x.ti+1 < truth_lens[i]) { // DEL
                        if (!done[x.hi][x.cri][x.ti+1] &&
                                !contains(done_this_wave,
                                    idx1(x.hi, x.cri, x.ti+1))) {
                            queue.push(idx1(x.hi, x.cri, x.ti+1));
                            done_this_wave.insert(idx1(x.hi, x.cri, x.ti+1));
                        }
                        if (!done[x.hi][x.cri][x.ti+1])
                            ptrs[x.hi][x.cri][x.ti+1] |= PTR_LEFT;
                    }
                    if (x.cri+1 < ref_len && x.ti+1 < truth_lens[i]) { // SUB
                        if (!done[x.hi][x.cri+1][x.ti+1] &&
                                !contains(done_this_wave, 
                                    idx1(x.hi, x.cri+1, x.ti+1))) {
                            queue.push(idx1(x.hi, x.cri+1, x.ti+1));
                            done_this_wave.insert(idx1(x.hi, x.cri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.cri+1][x.ti+1])
                            ptrs[x.hi][x.cri+1][x.ti+1] |= PTR_SUB;
                    }
                }
            }
            this_wave.clear();
            s[i]++;
        } // while loop (this alignment)

        // save where to start backtrack (prefer ref: omit vars which don't reduce ED)
        if (done[ri][ref_len-1][truth_lens[i]-1]) {
            pr_query_ref_end[i] = ri;
        } else if (done[ci][query_lens[i]-1][truth_lens[i]-1]) {
            pr_query_ref_end[i] = ci;
        } else { ERROR("Alignment not finished in 'prec_recall_aln()'."); }

    } // 4 alignments
} // function


/******************************************************************************/


int store_phase( 
        std::shared_ptr<superclusterData> clusterdata_ptr, 
        std::string ctg,
        std::vector<int> & s
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
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, std::string ctg,
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_scores, 
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs,
        std::vector<int> pr_query_ref_end
        ) {

    std::vector< std::vector<int> > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector<int> > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector< std::vector<int> > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector<int> pr_query_ref_beg(4);
    std::vector< std::vector<bool> > ref_loc_sync;

    for (int i = 0; i < 4; i++) {

        // init
        int ci = i*2 + QUERY;
        int ri = i*2 + REF;
        path.push_back(std::vector<idx1>());
        sync.push_back(std::vector<bool>());
        edits.push_back(std::vector<bool>());
        ref_loc_sync.push_back(std::vector<bool>(int(aln_ptrs[ri].size()), true));
        path_ptrs.push_back(std::vector< std::vector<int> >(aln_ptrs[ci].size(), 
            std::vector<int>(aln_ptrs[ci][0].size(), PTR_NONE)));
        path_ptrs.push_back(std::vector< std::vector<int> >(aln_ptrs[ri].size(), 
            std::vector<int>(aln_ptrs[ri][0].size(), PTR_NONE)));
        path_scores.push_back(std::vector< std::vector<int> >(aln_ptrs[ci].size(), 
            std::vector<int>(aln_ptrs[ci][0].size(), -1)));
        path_scores.push_back(std::vector< std::vector<int> >(aln_ptrs[ri].size(), 
            std::vector<int>(aln_ptrs[ri][0].size(), -1)));

        // backtrack start
        std::queue<idx1> queue;
        int start_hi = pr_query_ref_end[i];
        int start_cri = aln_ptrs[start_hi].size()-1;
        int start_ti = aln_ptrs[start_hi][0].size()-1;
        idx1 start(start_hi, start_cri, start_ti);
        path_ptrs[start_hi][start_cri][start_ti] = PTR_DIAG;
        aln_ptrs[start_hi][start_cri][start_ti] |= MAIN_PATH;
        path_scores[start_hi][start_cri][start_ti] = 0;
        queue.push(start);

        while (!queue.empty()) {
            idx1 x = queue.front(); queue.pop(); // current cell

            // ref locations which consume a query base and aren't 
            // on the main diagonal cannot be sync points
            if (aln_ptrs[x.hi][x.cri][x.ti] & (PTR_DIAG | PTR_SUB | PTR_INS)) {
                if (x.hi == ri && x.cri != truth_ref_ptrs[i][x.ti])
                    ref_loc_sync[i][x.cri] = false;
                if (x.hi == ci && query_ref_ptrs[i][x.cri] >= 0 && 
                        query_ref_ptrs[i][x.cri] != truth_ref_ptrs[i][x.ti])
                    ref_loc_sync[i][query_ref_ptrs[i][x.cri]] = false;
            }

            if (aln_ptrs[x.hi][x.cri][x.ti] & PTR_DIAG && // MATCH
                    x.cri > 0 && x.ti > 0) {

                // add to path
                idx1 y = idx1(x.hi, x.cri-1, x.ti-1); // next cell
                aln_ptrs[y.hi][y.cri][y.ti] |= PATH;

                // check for fp
                int is_fp = FALSE;
                if (x.hi == ri && ref_query_ptrs[i][y.cri] >= 0 &&
                        ref_query_ptrs[i][x.cri] != ref_query_ptrs[i][y.cri]+1) { // fp
                    is_fp = TRUE;
                }

                // update score
                if (path_scores[x.hi][x.cri][x.ti]+is_fp > 
                        path_scores[y.hi][y.cri][y.ti]) {
                    if (y.cri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.cri][y.ti] = PTR_DIAG;
                    path_scores[y.hi][y.cri][y.ti] = 
                        path_scores[x.hi][x.cri][x.ti]+is_fp;
                    queue.push(y);
                }

            }

            if (aln_ptrs[x.hi][x.cri][x.ti] & PTR_SUB && // SUB
                    x.cri > 0 && x.ti > 0) {

                // add to path
                idx1 y = idx1(x.hi, x.cri-1, x.ti-1);
                aln_ptrs[y.hi][y.cri][y.ti] |= PATH;

                // check for fp
                int is_fp = FALSE;
                if (x.hi == ri && ref_query_ptrs[i][y.cri] >= 0 &&
                        ref_query_ptrs[i][x.cri] != ref_query_ptrs[i][y.cri]+1) { // fp
                    is_fp = TRUE;
                }

                // update score
                if (path_scores[x.hi][x.cri][x.ti]+is_fp > 
                        path_scores[y.hi][y.cri][y.ti]) {
                    if (y.cri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.cri][y.ti] = PTR_SUB;
                    path_scores[y.hi][y.cri][y.ti] = 
                        path_scores[x.hi][x.cri][x.ti]+is_fp;
                    queue.push(y);
                }
            }

            if (aln_ptrs[x.hi][x.cri][x.ti] & PTR_INS && x.cri > 0) { // INS

                // add to path
                idx1 y = idx1(x.hi, x.cri-1, x.ti);
                aln_ptrs[y.hi][y.cri][y.ti] |= PATH;

                // check for fp
                int is_fp = FALSE;
                if (x.hi == ri && ref_query_ptrs[i][y.cri] >= 0 &&
                        ref_query_ptrs[i][x.cri] != ref_query_ptrs[i][y.cri]+1) { // fp
                    is_fp = TRUE;
                }

                // update score
                if (path_scores[x.hi][x.cri][x.ti]+is_fp > 
                        path_scores[y.hi][y.cri][y.ti]) {
                    if (y.cri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.cri][y.ti] = PTR_INS;
                    path_scores[y.hi][y.cri][y.ti] = 
                        path_scores[x.hi][x.cri][x.ti]+is_fp;
                    queue.push(y);
                } 
            }

            if (aln_ptrs[x.hi][x.cri][x.ti] & PTR_DEL && x.ti > 0) { // DEL

                // add to path
                idx1 y = idx1(x.hi, x.cri, x.ti-1);
                aln_ptrs[y.hi][y.cri][y.ti] |= PATH;

                // update score
                if (path_scores[x.hi][x.cri][x.ti] > 
                        path_scores[y.hi][y.cri][y.ti]) {
                    if (y.cri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                    path_ptrs  [y.hi][y.cri][y.ti] = PTR_DEL;
                    path_scores[y.hi][y.cri][y.ti] = 
                        path_scores[x.hi][x.cri][x.ti];
                    queue.push(y);
                }
            }

            if (x.hi == ri) { // REF -> QUERY SWAP
                if (aln_ptrs[x.hi][x.cri][x.ti] & PTR_SWAP &&
                        ref_query_ptrs[i][x.cri] >= 0) {

                    // add to path
                    idx1 y = idx1(ci, ref_query_ptrs[i][x.cri], x.ti);
                    aln_ptrs[y.hi][y.cri][y.ti] |= PATH;

                    // update score
                    if (path_scores[x.hi][x.cri][x.ti] > 
                            path_scores[y.hi][y.cri][y.ti]) {
                        if (y.cri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                        path_ptrs  [y.hi][y.cri][y.ti] = PTR_SWAP;
                        path_scores[y.hi][y.cri][y.ti] = 
                            path_scores[x.hi][x.cri][x.ti];
                        queue.push(y);
                    }
                }

            } else { // QUERY -> REF SWAP
                if (aln_ptrs[x.hi][x.cri][x.ti] & PTR_SWAP &&
                        query_ref_ptrs[i][x.cri] >= 0) {

                    // add to path
                    idx1 y = idx1(ri, query_ref_ptrs[i][x.cri], x.ti);
                    aln_ptrs[y.hi][y.cri][y.ti] |= PATH;

                    // update score
                    if (path_scores[x.hi][x.cri][x.ti] > 
                            path_scores[y.hi][y.cri][y.ti]) {
                        if (y.cri == 0 && y.ti == 0) pr_query_ref_beg[i] = y.hi;
                        path_ptrs  [y.hi][y.cri][y.ti] = PTR_SWAP;
                        path_scores[y.hi][y.cri][y.ti] = 
                            path_scores[x.hi][x.cri][x.ti];
                        queue.push(y);
                    }
                }
            }
        }
    }

    // get path and sync points
    get_prec_recall_path_sync(path, sync, edits, ref_loc_sync,
            path_ptrs, aln_ptrs, pr_query_ref_beg,
            query1_ref_ptrs, ref_query1_ptrs, query2_ref_ptrs, ref_query2_ptrs,
            truth1_ref_ptrs, truth2_ref_ptrs
    );

}


/******************************************************************************/


void get_prec_recall_path_sync(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector<bool> > ref_loc_sync, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector<int> & pr_query_ref_beg,
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs
        ) {

    // query <-> ref pointers
    std::vector< std::vector<int> > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector<int> > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector< std::vector<int> > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };

    for (int i = 0; i < 4; i++) {

        // init
        int ci = i*2 + QUERY;
        int ri = i*2 + REF;

        // path start
        int hi = pr_query_ref_beg[i];
        int cri = 0;
        int ti = 0;

        // first position is sync point
        path[i].push_back(idx1(hi, cri, ti));
        sync[i].push_back(true);
        aln_ptrs[hi][cri][ti] |= MAIN_PATH | PTR_SYNC;
        path_ptrs[hi][cri][ti] |= MAIN_PATH;

        // follow best-path pointers
        while (cri < int(path_ptrs[hi].size()) || ti < int(path_ptrs[hi][0].size())) {
            if (path_ptrs[hi][cri][ti] & PTR_DIAG) {
                cri++; ti++; edits[i].push_back(false);
            } else if (path_ptrs[hi][cri][ti] & PTR_SUB) {
                cri++; ti++; edits[i].push_back(true);
            } else if (path_ptrs[hi][cri][ti] & PTR_INS) {
                cri++; edits[i].push_back(true);
            } else if (path_ptrs[hi][cri][ti] & PTR_DEL) {
                ti++; edits[i].push_back(true);
            } else if (path_ptrs[hi][cri][ti] & PTR_SWAP) {
                if (hi == ri) {
                    hi = ci;
                    cri = ref_query_ptrs[i][cri];
                } else if(hi == ci) { // hi == ci
                    hi = ri;
                    cri = query_ref_ptrs[i][cri];
                } else {
                    ERROR("Unexpected hap index (%d)", hi);
                } 
                edits[i].push_back(false);
            } else {
                ERROR("No pointer for MAIN_PATH at (%d, %d, %d)", hi, cri, ti);
            }

            if (cri >= int(aln_ptrs[hi].size()) ||
                    ti >= int(aln_ptrs[hi][0].size())) break;

            // add point to path
            path[i].push_back(idx1(hi, cri, ti));
            aln_ptrs[hi][cri][ti] |= MAIN_PATH;
            path_ptrs[hi][cri][ti] |= MAIN_PATH;

            // determine if sync point
            bool is_sync;
            if (hi == ci) {
                if (query_ref_ptrs[i][cri] >= 0 &&  // on main diag
                    query_ref_ptrs[i][cri] == truth_ref_ptrs[i][ti]) {
                    is_sync = ref_loc_sync[i][ query_ref_ptrs[i][cri] ];
                } else {
                    is_sync = false;
                }
            } else { // hi == ri
                if (cri == truth_ref_ptrs[i][ti]) { // on main diag
                    is_sync = ref_loc_sync[i][cri];
                } else {
                    is_sync = false;
                }
            }
            if (is_sync) aln_ptrs[hi][cri][ti] |= PTR_SYNC;
            sync[i].push_back(is_sync);
        }
        // last position is sync
        sync[i][sync[i].size()-1] = true;
        idx1 last = path[i][path[i].size()-1];
        aln_ptrs[last.hi][last.cri][last.ti] |= PTR_SYNC;
    }
}


/******************************************************************************/


void calc_prec_recall(
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, std::string ctg,
        std::string query1, std::string query2, 
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs,
        std::vector< std::vector<idx1> > & path,
        std::vector< std::vector<bool> > & sync,
        std::vector< std::vector<bool> > & edits,
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector<int> pr_query_ref_end, int phase, int print
        ) {

    // set query/truth strings and pointers
    int beg = clusterdata_ptr->ctg_superclusters[ctg]->begs[sc_idx];
    std::vector<std::string> query = {query1, query1, query2, query2};
    std::vector<std::string> truth = {truth1, truth2, truth1, truth2};
    std::vector< std::vector<int> > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector<int> > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector< std::vector<int> > truth_ref_ptrs = { 
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
        int ci = i*2 + QUERY; // call index
        int qhi = i >> 1;      // query hap index
        int thi = i & 1;      // truth hap index

        // set variant ranges
        std::shared_ptr<ctgVariants> query_vars = 
                clusterdata_ptr->ctg_superclusters[ctg]->ctg_variants[QUERY][qhi];
        std::shared_ptr<ctgVariants> truth_vars = 
                clusterdata_ptr->ctg_superclusters[ctg]->ctg_variants[TRUTH][thi];
        int query_beg_idx = query_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[QUERY][qhi][sc_idx]];
        int query_end_idx = query_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[QUERY][qhi][sc_idx+1]];
        int truth_beg_idx = truth_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[TRUTH][thi][sc_idx]];
        int truth_end_idx = truth_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[TRUTH][thi][sc_idx+1]];

        // debug print
        if (print) printf("Alignment %s, aln_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(aln_ptrs[ci], query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(aln_ptrs[ri], ref, truth[i]);
        if (print) printf("Alignment %s, path_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(path_ptrs[ci], query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(path_ptrs[ri], ref, truth[i]);

        // init
        int hi = pr_query_ref_end[i];
        int cri = aln_ptrs[hi].size()-1;
        int ti = aln_ptrs[hi][0].size()-1;
        int new_ed = 0;
        int query_var_ptr = query_end_idx-1;
        int query_var_pos = query_vars->poss[query_var_ptr] - beg;
        int prev_query_var_ptr = query_var_ptr;
        int truth_var_ptr = truth_end_idx-1;
        int truth_var_pos = truth_vars->poss[truth_var_ptr] - beg;
        int prev_truth_var_ptr = truth_var_ptr;
        int pidx = path[i].size()-1;

        if (print) printf("\n%s:\n", aln_strs[i].data());
        while (pidx >= 0) {

            // set FP if necessary
            int cref_pos = (hi == ri) ? cri : query_ref_ptrs[i][cri];
            if (cref_pos < 0) cref_pos = query_var_pos+1;
            while (cref_pos < query_var_pos) {
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

            while (truth_ref_ptrs[i][ti] >= 0 && 
                    truth_ref_ptrs[i][ti] < truth_var_pos) { // passed REF variant
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
                        case TYPE_GRP:
                            old_ed += truth_vars->alts[truth_var_idx].length();
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
                        (hi == ri) ? "REF" : "QUERY", cri, ti,
                        (hi == ri) ? cri : query_ref_ptrs[i][cri],
                        truth_ref_ptrs[i][ti], old_ed, new_ed, 
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
                    (hi == ri) ? "REF  " : "QUERY", cri, ti,
                    edits[i][pidx] ? "X" : "=");

            // update path pointer
            pidx--;
            if (pidx < 0) break;

            // update location
            cri = path[i][pidx].cri;
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
        std::string query1, std::string query2, 
        std::string truth1, std::string truth2,
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
    INFO("Calculating edit distance");

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
            std::string query1 = "", ref_c1 = ""; 
            std::vector<int> query1_ref_ptrs, ref_query1_ptrs;
            generate_ptrs_strs(
                    query1, ref_c1,
                    query1_ref_ptrs, ref_query1_ptrs, 
                    sc->ctg_variants[QUERY][HAP1], sc->ctg_variants[QUERY][HAP1],
                    sc->superclusters[QUERY][HAP1][sc_idx], 0,
                    sc->superclusters[QUERY][HAP1][sc_idx+1], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string query2 = "", ref_c2 = ""; 
            std::vector<int> query2_ref_ptrs, ref_query2_ptrs;
            generate_ptrs_strs(
                    query2, ref_c2,
                    query2_ref_ptrs, ref_query2_ptrs, 
                    sc->ctg_variants[QUERY][HAP2], sc->ctg_variants[QUERY][HAP2],
                    sc->superclusters[QUERY][HAP2][sc_idx], 0,
                    sc->superclusters[QUERY][HAP2][sc_idx+1], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth1 = "", ref_t1 = ""; 
            std::vector<int> truth1_ref_ptrs, ref_truth1_ptrs;
            generate_ptrs_strs(
                    truth1, ref_t1,
                    truth1_ref_ptrs, ref_truth1_ptrs, 
                    sc->ctg_variants[TRUTH][HAP1], sc->ctg_variants[TRUTH][HAP1],
                    sc->superclusters[TRUTH][HAP1][sc_idx], 0,
                    sc->superclusters[TRUTH][HAP1][sc_idx+1], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth2 = "", ref_t2 = ""; 
            std::vector<int> truth2_ref_ptrs, ref_truth2_ptrs;
            generate_ptrs_strs(
                    truth2, ref_t2,
                    truth2_ref_ptrs, ref_truth2_ptrs, 
                    sc->ctg_variants[TRUTH][HAP2], sc->ctg_variants[TRUTH][HAP2],
                    sc->superclusters[TRUTH][HAP2][sc_idx], 0,
                    sc->superclusters[TRUTH][HAP2][sc_idx+1], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );

            // calculate four forward-pass alignment edit dists
            // query1-truth2, query1-truth1, query2-truth1, query2-truth2
            std::vector<int> aln_score(4), aln_query_ref_end(4);
            std::vector< std::vector< std::vector<int> > > aln_ptrs;
            calc_prec_recall_aln(
                    query1, query2, truth1, truth2, ref_c1,
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs,
                    aln_score, aln_ptrs, aln_query_ref_end
            );

            // store optimal phasing for each supercluster
            // ORIG: query1-truth1 and query2-truth2
            // SWAP: query1-truth2 and query2-truth1
            int phase = store_phase(clusterdata_ptr, ctg, aln_score);

            // calculate paths from alignment
            std::vector< std::vector<idx1> > path;
            std::vector< std::vector<bool> > sync;
            std::vector< std::vector<bool> > edit;
            std::vector< std::vector< std::vector<int> > > path_ptrs;
            std::vector< std::vector< std::vector<int> > > path_scores;
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
                    clusterdata_ptr, sc_idx, ctg,
                    query1, query2, truth1, truth2, ref_c1,
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs,
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    path, sync, edit, aln_ptrs, path_ptrs, 
                    aln_query_ref_end, phase, g.print_verbosity >= 2
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

            // calculate quality thresholds 
            // (where string would change when including/excluding variants)
            std::set<int> quals = {};
            for (int hap = 0; hap < HAPS; hap++) {
                int beg_idx = sc->ctg_variants[QUERY][hap]->clusters[
                        sc->superclusters[QUERY][hap][sc_idx]];
                int end_idx = sc->ctg_variants[QUERY][hap]->clusters[
                        sc->superclusters[QUERY][hap][sc_idx+1]];
                for (int var_idx = beg_idx; var_idx < end_idx; var_idx++) {
                    quals.insert(sc->ctg_variants[QUERY][hap]->var_quals[var_idx]+1);
                }
            }
            quals.insert(g.max_qual+2);

            // phasing is known, add scores for each hap
            std::vector<int> qual_dists(g.max_qual+2, 0);
            for (int hap = 0; hap < HAPS; hap++) {

                // sweep through quality thresholds
                int prev_qual = 0;
                for (int qual : quals) {

                    // generate query string (only applying variants with Q>=qual)
                    std::string query = generate_str(
                            clusterdata_ptr->ref, 
                            sc->ctg_variants[QUERY][hap], 
                            ctg, 
                            sc->ctg_variants[QUERY][hap]->clusters[
                                sc->superclusters[QUERY][hap][sc_idx]],
                            sc->ctg_variants[QUERY][hap]->clusters[
                                sc->superclusters[QUERY][hap][sc_idx+1]],
                            sc->begs[sc_idx], 
                            sc->ends[sc_idx], 
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
            if (g.print_verbosity >= 1) {
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

                printf("ORIG:      %s\n", ref_c1.data());
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


/******************************************************************************/


/* Perform Djikstra Smith-Waterman alignment of two strings, returning 
 * the farthest-reaching reference index of lesser or equal score to that provided.
 * Strings are generated by applying variants to draft ref, and skipping 
 * variants is not allowed. Neither is a diagonal ref transition.
 */
int sw_max_reach(std::string query, std::string ref, 
        std::vector<int> query_ref_ptrs,
        std::vector<int> ref_query_ptrs,
        int sub, int open, int extend,
        int score, bool reverse /*= false*/) {
    
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
        /* printf("s = %d\n", s); */

        // EXTEND WAVEFRONT (stay at same score)
        while (!queue.empty()) {
            idx2 x = queue.front(); queue.pop();
            /* printf("  x = (%c, %d, %d)\n", std::string("SID")[x.mi], x.ci, x.ri); */
            waves[s].insert(x);

            // allow non-diagonal match
            if (x.ci+1 < query_len && x.ri+1 < ref_len && ! (
                        query_ref_ptrs[x.ci] == x.ri &&
                        ref_query_ptrs[x.ri] == x.ci) &&
                    query[x.ci+1] == ref[x.ri+1] && 
                    !contains(done, {x.mi, x.ci+1, x.ri+1})) {
                idx2 next(x.mi, x.ci+1, x.ri+1);
                queue.push(next);
                done.insert(next);
            }

            // allow exiting D/I state (no penalty for forward only)
            if (!reverse && (x.mi == MAT_INS || x.mi == MAT_DEL) && 
                    !contains(done, {MAT_SUB, x.ci, x.ri})) {
                queue.push({MAT_SUB, x.ci, x.ri});
                done.insert({MAT_SUB, x.ci, x.ri});
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
                if (x.mi == MAT_SUB && x.ci+1 < query_len && x.ri+1 < ref_len) {
                    idx2 next({x.mi, x.ci+1, x.ri+1});
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
                if (x.mi == MAT_SUB && x.ci+1 < query_len) {
                    idx2 next({MAT_INS, x.ci+1, x.ri});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }

                // DEL opening
                if (x.mi == MAT_SUB && x.ri+1 < ref_len) {
                    idx2 next({MAT_DEL, x.ci, x.ri+1});
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
                        !contains(done, {MAT_SUB, x.ci, x.ri})) {
                    queue.push({MAT_SUB, x.ci, x.ri});
                    done.insert({MAT_SUB, x.ci, x.ri});
                }
            }
        }

        // INS/DEL extension
        extend_wave = s + 1 - extend;
        if (extend_wave >= 0 && extend_wave <= score) {
            for(idx2 x : waves[extend_wave]) {

                // INS extending (only from MAT_INS)
                if (x.mi == MAT_INS && x.ci+1 < query_len) {
                    idx2 next({x.mi, x.ci+1, x.ri});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }

                // DEL extending (only from MAT_DEL)
                if (x.mi == MAT_DEL && x.ri+1 < ref_len) {
                    idx2 next({x.mi, x.ci, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }
            }
        }
    } // end reach

    // search for farthest reach
    int min_wave = std::max(0, std::min(std::min(
                    open_wave, extend_wave), sub_wave));
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
            for (size_t cluster = 0; cluster < vars->clusters.size()-1; cluster++) {
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
                std::unordered_map<idx2, idx2> ptrs = sw_align(query, ref,
                        sub, open, extend);
                
                // backtrack
                std::vector<int> cigar = sw_backtrack(query, ref, ptrs);
                
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
    std::vector< std::unordered_set<idx2> > waves;

    // perform global alignment
    int s = 0;
    while(true) {
        waves.push_back(std::unordered_set<idx2>());
        /* printf("s = %d\n", s); */

        // EXTEND WAVEFRONT (stay at same score)
        while (!queue.empty()) {
            idx2 x = queue.front(); queue.pop();
            /* printf("  x = (%c, %d, %d)\n", std::string("SID")[x.mi], x.ci, x.ri); */
            waves[s].insert(x);

            // allow non-diagonal match
            if (x.mi == MAT_SUB &&
                    x.ci+1 < query_len && x.ri+1 < ref_len &&
                    query[x.ci+1] == ref[x.ri+1] && 
                    !contains(done, {x.mi, x.ci+1, x.ri+1})) {
                idx2 next(x.mi, x.ci+1, x.ri+1);
                queue.push(next); done.insert(next); ptrs[next] = x;
            }

            // allow exiting D/I state freely
            if ((x.mi == MAT_INS || x.mi == MAT_DEL) && 
                    !contains(done, {MAT_SUB, x.ci, x.ri})) {
                idx2 next(MAT_SUB, x.ci, x.ri);
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
                if (x.mi == MAT_SUB && x.ci+1 < query_len && x.ri+1 < ref_len) {
                    idx2 next(x.mi, x.ci+1, x.ri+1);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;
                    }
                }
            }
        }

        // INS/DEL opening (only from MAT_SUB)
        open_wave = s + 1 - (open+extend);
        if (open_wave >= 0) {
            for(idx2 x : waves[open_wave]) {

                // INS opening
                if (x.mi == MAT_SUB && x.ci+1 < query_len) {
                    idx2 next(MAT_INS, x.ci+1, x.ri);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;
                    }
                }

                // DEL opening
                if (x.mi == MAT_SUB && x.ri+1 < ref_len) {
                    idx2 next(MAT_DEL, x.ci, x.ri+1);
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;
                    }
                }
            }
        }

        // INS/DEL extension
        extend_wave = s + 1 - extend;
        if (extend_wave >= 0) {
            for(idx2 x : waves[extend_wave]) {

                // INS extending (only from MAT_INS)
                if (x.mi == MAT_INS && x.ci+1 < query_len) {
                    idx2 next({x.mi, x.ci+1, x.ri});
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;
                    }
                }

                // DEL extending (only from MAT_DEL)
                if (x.mi == MAT_DEL && x.ri+1 < ref_len) {
                    idx2 next({x.mi, x.ci, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next); done.insert(next); ptrs[next] = x;
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

    while (pos.ci >= 0 || pos.ri >= 0) {
        auto next_itr = ptrs.find(pos);
        if (next_itr == ptrs.end())
            ERROR("(%d, %d, %d) not found in ptrs during backtrack.",
                    pos.mi, pos.ci, pos.ri);
        idx2 next = ptrs.find(pos)->second;
        if (next.mi == pos.mi) { // stay in same matrix
            switch (pos.mi) {
                case MAT_SUB:
                    if (query[pos.ci] == ref[pos.ri]) { // match
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
            }
        }
        pos = next;
    }
    return cigar;
}
