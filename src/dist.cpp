#include <string>
#include <vector>
#include <cstdio>
#include <chrono>
#include <utility>
#include <queue>
#include <unordered_set>

#include "dist.h"
#include "print.h"
#include "cluster.h"


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
                                    TYPE_SUB, INSIDE, std::string(1,ref_str[ref_idx]), 
                                    std::string(1,alt_str[alt_idx]), GT_REF_REF, 60, 60);
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
                                    indel_len, h, TYPE_DEL, INSIDE,
                                    ref_str.substr(ref_idx, indel_len),
                                    "", GT_REF_REF, 60, 60);
                            ref_idx += indel_len;
                            break;

                        case PTR_UP: // insertion
                            cig_idx++; indel_len++;

                            // multi-base insertion
                            while (cig_idx < cig.size() && cig[cig_idx] == PTR_UP) {
                                cig_idx++; indel_len++;
                            }
                            results.ctg_variants[h][ctg]->add_var(beg+ref_idx,
                                    0, h, TYPE_INS, INSIDE, "", 
                                    alt_str.substr(alt_idx, indel_len), GT_REF_REF, 60, 60);
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


/* For each position on hap1 and hap2, it generates pointers from one hap to the
 * other if both match the reference, or else -1. This method also generates the
 * strings for alignment, and the pretty-printing strings for debug output.
 */
void generate_ptrs_strs(
        std::string & hap1, std::string & hap2,          // actual strings 
        std::string & hap1_str, std::string & hap2_str,  // colored debug strs
        std::vector<int> & hap1_ptrs, std::vector<int> & hap2_ptrs,
        std::shared_ptr<ctgVariants> hap1_vars, std::shared_ptr<ctgVariants> hap2_vars,
        size_t hap1_clust_beg_idx, size_t hap2_clust_beg_idx,
        size_t hap1_clust_end_idx, size_t hap2_clust_end_idx,
        int beg_pos, int end_pos, std::shared_ptr<fastaData> ref, std::string ctg
        ) {

    // generate hap1 and hap2 strings and pointers
    int hap1_var_idx = hap1_vars->clusters[hap1_clust_beg_idx];
    int hap2_var_idx = hap2_vars->clusters[hap2_clust_beg_idx];
    for (int hap1_ref_pos = beg_pos, hap2_ref_pos = beg_pos; 
            hap1_ref_pos < end_pos || hap2_ref_pos < end_pos; ) {

        // CONSIDER TRUTH1 ONLY, PRIOR REFERENCE POSITION
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
                    ERROR("Contig %s not present in reference FASTA",
                            ctg.data());
                    exit(1);
                }
            }
        }

        // CONSIDER TRUTH2 ONLY, PRIOR REFERENCE POSITION
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
                    ERROR("Contig %s not present in reference FASTA",
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

            // BOTH MATCH REFERENCE, ADD POINTERS FOR TRANSITIONS
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
                    ERROR("Contig '%s' not present in reference FASTA",
                            ctg.data());
                    exit(1);
                }
            }
        }
    }
}


/******************************************************************************/


class idx {
public:
    int hi;  // hap idx
    int cri; // calls/ref idx
    int ti;  // truth idx

    idx() : hi(0), cri(0), ti(0) {};
    idx(int h, int c, int t) : hi(h), cri(c), ti(t) {};
    bool operator<(const idx & idx2) const {
        if (this->hi < idx2.hi) return true;
        else if (this->hi == idx2.hi && this->cri < idx2.cri) return true;
        else if (this->hi == idx2.hi && this->cri == idx2.cri && this->ti < idx2.ti) return true;
        return false;
    }
    bool operator==(const idx & idx2) const {
        return this->hi == idx2.hi && this->cri == idx2.cri && this->ti == idx2.ti;
    }
};

namespace std {
    template<> struct hash<idx> {
        std::uint64_t operator()(const idx& idx1) const noexcept {
            return uint64_t(idx1.hi)<<60 | uint64_t(idx1.cri) << 30 | uint64_t(idx1.ti);
        }
    };
}

bool contains(const std::unordered_set<idx> & wave, const idx & idx) {
    return wave.find(idx) != wave.end();
}


void calc_prec_recall_aln(
        std::string calls1, std::string calls2,
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> calls1_ref_ptrs, std::vector<int> ref_calls1_ptrs,
        std::vector<int> calls2_ref_ptrs, std::vector<int> ref_calls2_ptrs,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & ptrs,
        std::vector<int> & pr_calls_ref_end
        ) {
    
    // set loop variables
    int ref_len = ref.size();
    std::vector<std::string> calls {calls1, calls1, calls2, calls2};
    std::vector<std::string> truth {truth1, truth2, truth1, truth2};
    std::vector< std::vector<int> > calls_ref_ptrs {
            calls1_ref_ptrs, calls1_ref_ptrs, calls2_ref_ptrs, calls2_ref_ptrs };
    std::vector< std::vector<int> > ref_calls_ptrs {
            ref_calls1_ptrs, ref_calls1_ptrs, ref_calls2_ptrs, ref_calls2_ptrs };
    std::vector<int> calls_lens = 
            {int(calls1.size()), int(calls1.size()), int(calls2.size()), int(calls2.size())};
    std::vector<int> truth_lens = 
            {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};

    std::vector< std::vector< std::vector<bool> > > done;

    // for each combination of calls and truth
    for (int i = 0; i < 4; i++) {
        int ci = 2*i + CALLS; // calls index (offs and ptrs)
        int ri = 2*i + REF;   // ref index   (offs and ptrs)

        // init full pointer/done matrices
        ptrs.push_back(std::vector< std::vector<int> >(calls_lens[i], 
                std::vector<int>(truth_lens[i], PTR_NONE)));
        ptrs.push_back(std::vector< std::vector<int> >(ref_len, 
                std::vector<int>(truth_lens[i], PTR_NONE)));
        ptrs[ci][0][0] = PTR_DIAG;
        ptrs[ri][0][0] = PTR_DIAG;

        done.push_back(std::vector< std::vector<bool> >(calls_lens[i], 
                std::vector<bool>(truth_lens[i], false)));
        done.push_back(std::vector< std::vector<bool> >(ref_len, 
                std::vector<bool>(truth_lens[i], false)));
        
        // set first wavefront
        std::queue< idx > queue;
        queue.push({ci, 0, 0});
        queue.push({ri, 0, 0});
        ptrs[ri][0][0] |= PTR_DIAG;
        done[ri][0][0] = true;
        ptrs[ci][0][0] |= PTR_DIAG;
        done[ci][0][0] = true;

        // continue looping until full alignment found
        while (!done[ci][calls_lens[i]-1][truth_lens[i]-1] &&
                !done[ri][ref_len-1][truth_lens[i]-1]) {
            /* printf("s = %d\n", s[i]); */
            if (queue.empty()) ERROR("Empty queue in 'prec_recall_aln()'.");

            // EXTEND WAVEFRONT (stay at same score)
            std::unordered_set< idx > this_wave;
            std::unordered_set< idx > done_this_wave;
            while (!queue.empty()) {
                idx x = queue.front(); queue.pop();                  // next cell
                /* printf("x = (%d, %d, %d)\n", x.hi, x.cri, x.ti); */
                this_wave.insert(x);                                       // record visit
                if (x.hi == ci) { // == CALLS
                    // allow match
                    if (calls[i][x.cri] == truth[i][x.ti] &&              // if match
                            x.cri+1 < calls_lens[i] &&                    // if in range
                            x.ti+1  < truth_lens[i]) {
                        if (!done[x.hi][x.cri+1][x.ti+1] && 
                                !contains(done_this_wave, idx(x.hi, x.cri+1,x.ti+1))) {
                            queue.push(idx(x.hi, x.cri+1, x.ti+1));
                            done_this_wave.insert(idx(x.hi, x.cri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.cri+1][x.ti+1])
                            ptrs[x.hi][x.cri+1][x.ti+1] |= PTR_DIAG;
                    }
                    // allow phase swap
                    if (calls_ref_ptrs[i][x.cri] >= 0) {                  // if swap allowed
                        if (!done[ri][calls_ref_ptrs[i][x.cri]][x.ti] &&
                                !contains(done_this_wave, idx(ri, calls_ref_ptrs[i][x.cri], x.ti))) {
                            queue.push(idx(ri, 
                                    calls_ref_ptrs[i][x.cri], 
                                    x.ti));
                            done_this_wave.insert(idx(ri, calls_ref_ptrs[i][x.cri], x.ti));
                        }
                        if (!done[ri][calls_ref_ptrs[i][x.cri]][x.ti])
                            ptrs[ri][calls_ref_ptrs[i][x.cri]][x.ti] |= PTR_SWAP;
                    }
                } else { // x.hi == ri == REF
                    // allow match
                    if (ref[x.cri] == truth[i][x.ti] &&                     // if match
                            x.cri+1 < ref_len &&                            // if in range
                            x.ti+1  < truth_lens[i]) {
                        if (!done[x.hi][x.cri+1][x.ti+1] &&
                                !contains(done_this_wave, idx(x.hi, x.cri+1, x.ti+1))) {
                            queue.push(idx(x.hi, x.cri+1, x.ti+1));
                            done_this_wave.insert(idx(x.hi, x.cri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.cri+1][x.ti+1])
                            ptrs[x.hi][x.cri+1][x.ti+1] |= PTR_DIAG;
                    }
                    // allow phase swap
                    if (ref_calls_ptrs[i][x.cri] >= 0) {                    // if swap allowed
                        if (!done[ci][ref_calls_ptrs[i][x.cri]][x.ti] &&
                                !contains(this_wave, idx(ci, ref_calls_ptrs[i][x.cri], x.ti))) {
                            queue.push(idx(ci, 
                                        ref_calls_ptrs[i][x.cri], 
                                        x.ti));
                            done_this_wave.insert(idx(ci, ref_calls_ptrs[i][x.cri], x.ti));
                        }
                        if (!done[ci][ref_calls_ptrs[i][x.cri]][x.ti])
                            ptrs[ci][ref_calls_ptrs[i][x.cri]][x.ti] |= PTR_SWAP;
                    }
                }
            }

            // mark all cells visited this wave as done
            for (auto x : done_this_wave) { done[x.hi][x.cri][x.ti] = true; }
            done_this_wave.clear();

            // exit if we're done aligning
            if (done[ci][calls_lens[i]-1][truth_lens[i]-1] ||
                done[ri][ref_len-1][truth_lens[i]-1]) break;


            // NEXT WAVEFRONT (increase score by one)
            for (auto x : this_wave) {
                if (x.hi == ci) {
                    if (x.cri+1 < calls_lens[i]) {
                        if (!done[x.hi][x.cri+1][x.ti] && 
                                !contains(done_this_wave,
                                    idx(x.hi, x.cri+1, x.ti))) {
                            queue.push(idx(x.hi, x.cri+1, x.ti));
                        }
                        if (!done[x.hi][x.cri+1][x.ti])
                            ptrs[x.hi][x.cri+1][x.ti] |= PTR_UP;
                    }
                    if (x.ti+1 < truth_lens[i]) {
                        if (!done[x.hi][x.cri][x.ti+1] &&
                                !contains(done_this_wave,
                                    idx(x.hi, x.cri, x.ti+1))) {
                            queue.push(idx(x.hi, x.cri, x.ti+1));
                        }
                        if (!done[x.hi][x.cri][x.ti+1])
                            ptrs[x.hi][x.cri][x.ti+1] |= PTR_LEFT;
                    }
                    if (x.cri+1 < calls_lens[i] && x.ti+1 < truth_lens[i]) {
                        if (!done[x.hi][x.cri+1][x.ti+1] &&
                                !contains(done_this_wave,
                                    idx(x.hi, x.cri+1, x.ti+1))) {
                            queue.push(idx(x.hi, x.cri+1, x.ti+1));
                        }
                        if (!done[x.hi][x.cri+1][x.ti+1])
                            ptrs[x.hi][x.cri+1][x.ti+1] |= PTR_SUB;
                    }
                } else { //x.hi == REF
                    if (x.cri+1 < ref_len) {
                        if (!done[x.hi][x.cri+1][x.ti] &&
                                !contains(done_this_wave,
                                    idx(x.hi, x.cri+1, x.ti))) {
                            queue.push(idx(x.hi, x.cri+1, x.ti));
                        }
                        if (!done[x.hi][x.cri+1][x.ti])
                            ptrs[x.hi][x.cri+1][x.ti] |= PTR_UP;
                    }
                    if (x.ti+1 < truth_lens[i]) {
                        if (!done[x.hi][x.cri][x.ti+1] &&
                                !contains(done_this_wave,
                                    idx(x.hi, x.cri, x.ti+1))) {
                            queue.push(idx(x.hi, x.cri, x.ti+1));
                        }
                        if (!done[x.hi][x.cri][x.ti+1])
                            ptrs[x.hi][x.cri][x.ti+1] |= PTR_LEFT;
                    }
                    if (x.cri+1 < ref_len && x.ti+1 < truth_lens[i]) {
                        if (!done[x.hi][x.cri+1][x.ti+1] &&
                                !contains(done_this_wave, 
                                    idx(x.hi, x.cri+1, x.ti+1))) {
                            queue.push(idx(x.hi, x.cri+1, x.ti+1));
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
            pr_calls_ref_end[i] = REF;
        } else if (done[ci][calls_lens[i]-1][truth_lens[i]-1]) {
            pr_calls_ref_end[i] = CALLS;
        } else { ERROR("Alignment not finished in 'prec_recall_aln()'."); }

    } // 4 alignments
} // function


/******************************************************************************/


int store_phase( 
        std::shared_ptr<clusterData> clusterdata_ptr, 
        std::string ctg,
        std::vector<int> & s
        ) {

    // calculate best phasing
    int orig_phase_dist = s[CALLS1_TRUTH1] + s[CALLS2_TRUTH2];
    int swap_phase_dist = s[CALLS2_TRUTH1] + s[CALLS1_TRUTH2];
    int phase = PHASE_NONE; // default either way if equal dist
    if (orig_phase_dist < swap_phase_dist) phase = PHASE_ORIG;
    if (swap_phase_dist < orig_phase_dist) phase = PHASE_SWAP;

    // save alignment information
    clusterdata_ptr->ctg_superclusters[ctg]->add_phasing(
            phase, orig_phase_dist, swap_phase_dist);
    return phase;
}


/******************************************************************************/


void color_ptrs(
        std::vector< std::vector< std::vector<int> > > & ptrs, 
        std::vector<int> calls1_ref_ptrs, std::vector<int> ref_calls1_ptrs,
        std::vector<int> calls2_ref_ptrs, std::vector<int> ref_calls2_ptrs,
        std::vector<int> pr_calls_ref_end
        ) {

    // calls <-> ref pointers
    std::vector< std::vector<int> > calls_ref_ptrs = { 
            calls1_ref_ptrs, calls1_ref_ptrs, calls2_ref_ptrs, calls2_ref_ptrs };
    std::vector< std::vector<int> > ref_calls_ptrs = { 
            ref_calls1_ptrs, ref_calls1_ptrs, ref_calls2_ptrs, ref_calls2_ptrs };

    for (int aln = 0; aln < 4; aln++) {

        // init
        int aln_idx = aln*2;
        int hap_idx = pr_calls_ref_end[aln];
        int prev_hap_idx = hap_idx;
        int ptr_refcalls = ptrs[aln_idx+hap_idx].size()-1;
        int ptr_truth = ptrs[aln_idx+hap_idx][0].size()-1;

        // RIGHT PATH
        while (ptr_refcalls >= 0 || ptr_truth >= 0) {
            ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] |= PTR_RPATH;
            // prefer REF (omit vars which don't reduce ED) but prevent loops
            if (hap_idx == CALLS && prev_hap_idx == CALLS && 
                    ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_SWAP) {
                ptr_refcalls = calls_ref_ptrs[aln][ptr_refcalls];
                if (ptr_refcalls < 0) ERROR("Backtracking OOB");
                prev_hap_idx = CALLS; hap_idx = REF; continue;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_UP) {
                ptr_refcalls--;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_DIAG) {
                ptr_refcalls--; ptr_truth--;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_SUB) {
                ptr_refcalls--; ptr_truth--;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_LEFT) {
                ptr_truth--;
            } else if (hap_idx == REF) {
                if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_SWAP) {
                    ptr_refcalls = ref_calls_ptrs[aln][ptr_refcalls];
                    if (ptr_refcalls < 0) ERROR("Backtracking OOB");
                    prev_hap_idx = REF; hap_idx = CALLS; continue;
                }
            } else {
                ERROR("Unexpected alignment pointer (%d) in color_ptrs()", 
                        ptrs[aln_idx + hap_idx][ptr_refcalls][ptr_truth]);
            }
            prev_hap_idx = hap_idx;
        }

        // LEFT PATH
        hap_idx = pr_calls_ref_end[aln];
        prev_hap_idx = hap_idx;
        ptr_refcalls = ptrs[aln_idx+hap_idx].size()-1;
        ptr_truth = ptrs[aln_idx+hap_idx][0].size()-1;
        while (ptr_refcalls >= 0 || ptr_truth >= 0) {
            ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] |= PTR_LPATH;
            // prefer REF (omit vars which don't reduce ED) but prevent loops
            if (hap_idx == CALLS && prev_hap_idx == CALLS &&
                    ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_SWAP) {
                ptr_refcalls = calls_ref_ptrs[aln][ptr_refcalls];
                if (ptr_refcalls < 0) ERROR("Backtracking OOB");
                prev_hap_idx = CALLS; hap_idx = REF; continue;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_LEFT) {
                ptr_truth--;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_DIAG) {
                ptr_refcalls--; ptr_truth--;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_SUB) {
                ptr_refcalls--; ptr_truth--;
            } else if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_UP) {
                ptr_refcalls--;
            } else if (hap_idx == REF) {
                if (ptrs[aln_idx+hap_idx][ptr_refcalls][ptr_truth] & PTR_SWAP) {
                    ptr_refcalls = ref_calls_ptrs[aln][ptr_refcalls];
                    if (ptr_refcalls < 0) ERROR("Backtracking OOB");
                    prev_hap_idx = REF; hap_idx = CALLS; continue;
                }
            } else {
                ERROR("Unexpected alignment pointer (%d) in color_ptrs()", 
                        ptrs[aln_idx + hap_idx][ptr_refcalls][ptr_truth]);
            }
            prev_hap_idx = hap_idx;
        }

    }
}


/******************************************************************************/


void calc_prec_recall(
        std::shared_ptr<clusterData> clusterdata_ptr, int sc_idx, std::string ctg,
        std::string calls1, std::string calls2, 
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> calls1_ref_ptrs, std::vector<int> ref_calls1_ptrs,
        std::vector<int> calls2_ref_ptrs, std::vector<int> ref_calls2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs,
        std::vector< std::vector< std::vector<int> > > & ptrs, 
        std::vector<int> pr_calls_ref_end, int phase, int print
        ) {

    // set useful vectors for indexing based on current alignment
    int beg = clusterdata_ptr->ctg_superclusters[ctg]->begs[sc_idx];
    std::shared_ptr<ctgVariants> calls1_vars = 
            clusterdata_ptr->ctg_superclusters[ctg]->calls1_vars;
    std::shared_ptr<ctgVariants> calls2_vars = 
            clusterdata_ptr->ctg_superclusters[ctg]->calls2_vars;
    std::shared_ptr<ctgVariants> truth1_vars = 
            clusterdata_ptr->ctg_superclusters[ctg]->truth1_vars;
    std::shared_ptr<ctgVariants> truth2_vars = 
            clusterdata_ptr->ctg_superclusters[ctg]->truth2_vars;
    std::vector<std::string> calls = {calls1, calls1, calls2, calls2};
    std::vector<std::string> truth = {truth1, truth2, truth1, truth2};
    std::vector< std::vector<int> > calls_ref_ptrs = { 
            calls1_ref_ptrs, calls1_ref_ptrs, calls2_ref_ptrs, calls2_ref_ptrs };
    std::vector< std::vector<int> > ref_calls_ptrs = { 
            ref_calls1_ptrs, ref_calls1_ptrs, ref_calls2_ptrs, ref_calls2_ptrs };
    std::vector< std::vector<int> > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector< std::shared_ptr<ctgVariants> > calls_vars = {
            calls1_vars, calls1_vars, calls2_vars, calls2_vars };
    std::vector< std::shared_ptr<ctgVariants> > truth_vars = {
            truth1_vars, truth2_vars, truth1_vars, truth2_vars };
    int calls1_beg_idx = calls1_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->calls1_beg_idx[sc_idx] ];
    int calls2_beg_idx = calls2_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->calls2_beg_idx[sc_idx] ];
    std::vector<int> calls_beg_idx = {
            calls1_beg_idx, calls1_beg_idx, calls2_beg_idx, calls2_beg_idx };
    int calls1_end_idx = calls1_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->calls1_end_idx[sc_idx] ];
    int calls2_end_idx = calls2_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->calls2_end_idx[sc_idx] ];
    std::vector<int> calls_end_idx = {
            calls1_end_idx, calls1_end_idx, calls2_end_idx, calls2_end_idx };

    // indices into ptr/off matrices depend on decided phasing
    std::vector<int> calls_indices;
    if (phase == PHASE_SWAP) {
        calls_indices.push_back(CALLS1_TRUTH2);
        calls_indices.push_back(CALLS2_TRUTH1);
    } else if (phase == PHASE_ORIG || phase == PHASE_NONE) { // keep
        calls_indices.push_back(CALLS1_TRUTH1);
        calls_indices.push_back(CALLS2_TRUTH2);
    } else {
        ERROR("Unexpected phase (%d)", phase);
    }

    // for only the selected phasing
    if (print) printf("\n=======================================================================\n");
    for (int aln : calls_indices) {

        // init
        int aln_idx = aln*2;
        int hap_idx = pr_calls_ref_end[aln];
        int ptr_refcalls = ptrs[aln_idx+hap_idx].size()-1;
        int ptr_truth = ptrs[aln_idx+hap_idx][0].size()-1;
        int new_ed = 0;
        bool main_diag = true;
        int calls_var_ptr = calls_end_idx[aln]-1;
        int calls_var_pos = calls_vars[aln]->poss[calls_var_ptr] - beg;
        int prev_calls_var_ptr = calls_var_ptr;

        // debug print pointer matrices
        if (print) print_ptrs(ptrs[aln_idx+CALLS], calls[aln], truth[aln]);
        if (print) print_ptrs(ptrs[aln_idx+REF], ref, truth[aln]);

        if (print) printf("\n%s:\n", aln_strs[aln].data());
        while (ptr_refcalls >= 0 || ptr_truth >= 0) {

            // set FP if necessary
            int ref_pos = (hap_idx == REF) ? ptr_refcalls : 
                    calls_ref_ptrs[aln][ptr_refcalls];
            if (ref_pos < calls_var_pos) {
                if (hap_idx == REF) { // FP
                    calls_vars[aln]->errtypes[calls_var_ptr] = ERRTYPE_FP;
                    calls_vars[aln]->credit[calls_var_ptr] = 0;
                    if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                            calls_vars[aln]->refs[calls_var_ptr].data(),
                            calls_vars[aln]->alts[calls_var_ptr].data(), "FP", 0.0f);
                    calls_var_ptr--;
                } else { // passed variant
                    calls_var_ptr--;
                }
            }

            // sync point: set TP/PP
            if (main_diag) {

                // calculate old edit distance
                int old_ed = 0;
                for (int calls_var_idx = prev_calls_var_ptr; 
                        calls_var_idx > calls_var_ptr; calls_var_idx--) {
                    switch (calls_vars[aln]->types[calls_var_idx]) {
                        case TYPE_SUB:
                            old_ed += 1;
                            break;
                        case TYPE_INS:
                            old_ed += calls_vars[aln]->alts[calls_var_idx].length();
                            break;
                        case TYPE_DEL:
                            old_ed += calls_vars[aln]->refs[calls_var_idx].length();
                            break;
                        case TYPE_GRP:
                            old_ed += calls_vars[aln]->alts[calls_var_idx].length();
                            old_ed += calls_vars[aln]->refs[calls_var_idx].length();
                            break;
                        default:
                            ERROR("Unexpected variant type (%d) in calc_prec_recall().",
                                    calls_vars[aln]->types[calls_var_idx]);
                            break;
                    }
                }

                for (int calls_var_idx = prev_calls_var_ptr; 
                        calls_var_idx > calls_var_ptr; calls_var_idx--) {
                    // don't overwrite FPs
                    if (calls_vars[aln]->errtypes[calls_var_idx] == ERRTYPE_UN) {
                        if (new_ed == 0) { // TP
                            calls_vars[aln]->errtypes[calls_var_idx] = ERRTYPE_TP;
                            calls_vars[aln]->credit[calls_var_idx] = 1;
                            if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                    calls_vars[aln]->refs[calls_var_idx].data(),
                                    calls_vars[aln]->alts[calls_var_idx].data(), 
                                    "TP", 1.0f);
                        } else { // PP
                            calls_vars[aln]->errtypes[calls_var_idx] = ERRTYPE_PP;
                            calls_vars[aln]->credit[calls_var_idx] = 1;
                            if (print) printf("REF='%s'\tALT='%s'\t%s\t%f\n",
                                    calls_vars[aln]->refs[calls_var_idx].data(),
                                    calls_vars[aln]->alts[calls_var_idx].data(), 
                                    "PP", 1 - (float(new_ed)/old_ed));
                        }
                    }
                }

                if (print) printf("SYNC @%s(%d,%d)=REF(%d,%d), ED %d->%d, VARS %d-%d\n",
                        (hap_idx == REF) ? "REF" : "CALLS", ptr_refcalls, ptr_truth,
                        (hap_idx == REF) ? ptr_refcalls : ref_calls_ptrs[aln][ptr_refcalls],
                        truth_ref_ptrs[aln][ptr_truth], old_ed, new_ed, 
                        prev_calls_var_ptr, calls_var_ptr);

                prev_calls_var_ptr = calls_var_ptr;
                new_ed = 0;
            }

            // update pointers and edit distance
            if (print) printf("%s %s (%d,%d)\n", 
                    main_diag ? "*" : " ",
                    hap_idx ? "REF  " : "CALLS", 
                    ptr_refcalls, ptr_truth);

            int ptr_val = ptrs[aln_idx + hap_idx][ptr_refcalls][ptr_truth];
            if (ptr_val & PTR_UP) {
                if (print) printf("INS\n");
                ptr_refcalls--;
                new_ed++;
            } else if (ptr_val & PTR_LEFT) {
                if (print) printf("DEL\n");
                ptr_truth--;
                new_ed++;
            } else if (ptr_val & PTR_DIAG) {
                ptr_refcalls--;
                ptr_truth--;
            } else if (ptr_val & PTR_SUB) {
                if (print) printf("SUB\n");
                ptr_refcalls--;
                ptr_truth--;
                new_ed++;
            } else if (ptr_val & PTR_SWAP) {
                if (print) printf("SWAP\n");
                if (hap_idx == CALLS) {
                    ptr_refcalls = calls_ref_ptrs[aln][ptr_refcalls];
                    if (ptr_refcalls < 0) ERROR("Backtracking OOB");
                    hap_idx = REF;
                } else if (hap_idx == REF) {
                    ptr_refcalls = ref_calls_ptrs[aln][ptr_refcalls];
                    if (ptr_refcalls < 0) ERROR("Backtracking OOB");
                    hap_idx = CALLS;
                } else {
                    ERROR("Unexpected hap (%d)", hap_idx);
                }
            } else {
                ERROR("Unexpected alignment pointer (%d) in calc_prec_recall()", 
                        ptrs[aln_idx + hap_idx][ptr_refcalls][ptr_truth]);
            }

            // update if on main diag
            main_diag = (hap_idx == REF) ? 
                    (ptr_refcalls == truth_ref_ptrs[aln][ptr_truth]) : // REF
                    (calls_ref_ptrs[aln][ptr_refcalls] == 
                     truth_ref_ptrs[aln][ptr_truth] && 
                     calls_ref_ptrs[aln][ptr_refcalls] >= 0); // CALLS

            // update next variant position
            calls_var_pos = (calls_var_ptr < calls_beg_idx[aln]) ? -1 :
                calls_vars[aln]->poss[calls_var_ptr] - beg;
        }
    }
}


/******************************************************************************/


void calc_edit_dist_aln(
        std::string calls1, std::string calls2, 
        std::string truth1, std::string truth2,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & offs,
        std::vector< std::vector< std::vector<int> > > & ptrs
        ) {

    // ALIGNMENT
    std::vector<std::string> calls {calls1, calls1, calls2, calls2};
    std::vector<std::string> truth {truth1, truth2, truth1, truth2};
    std::vector<int> calls_lens = 
            {int(calls1.size()), int(calls1.size()), int(calls2.size()), int(calls2.size())};
    std::vector<int> truth_lens = 
            {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};

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
}


/******************************************************************************/


void alignment_wrapper(std::shared_ptr<clusterData> clusterdata_ptr) {

    int distance = 0;
    for (std::string ctg : clusterdata_ptr->contigs) {

        // set variant pointers
        std::shared_ptr<ctgVariants> calls1_vars = clusterdata_ptr->ctg_superclusters[ctg]->calls1_vars;
        std::shared_ptr<ctgVariants> calls2_vars = clusterdata_ptr->ctg_superclusters[ctg]->calls2_vars;
        std::shared_ptr<ctgVariants> truth1_vars = clusterdata_ptr->ctg_superclusters[ctg]->truth1_vars;
        std::shared_ptr<ctgVariants> truth2_vars = clusterdata_ptr->ctg_superclusters[ctg]->truth2_vars;

        // set superclusters pointer
        std::shared_ptr<ctgClusters> sc = clusterdata_ptr->ctg_superclusters[ctg];

        // iterate over superclusters
        for(int sc_idx = 0; sc_idx < clusterdata_ptr->ctg_superclusters[ctg]->n; sc_idx++) {

            // PRECISION-RECALL
            std::string calls1_c1 = "", ref_c1 = "", calls1_str_c1 = "", ref_str_c1 = ""; 
            std::vector<int> calls1_ref_ptrs, ref_calls1_ptrs;
            generate_ptrs_strs( // calls1_vars[0:0] contains no variants -> ref
                    calls1_c1, ref_c1, calls1_str_c1, ref_str_c1, 
                    calls1_ref_ptrs, ref_calls1_ptrs, calls1_vars, calls1_vars,
                    sc->calls1_beg_idx[sc_idx], 0,
                    sc->calls1_end_idx[sc_idx], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string calls2_c2 = "", ref_c2 = "", calls2_str_c2 = "", ref_str_c2 = ""; 
            std::vector<int> calls2_ref_ptrs, ref_calls2_ptrs;
            generate_ptrs_strs( // calls2_vars[0:0] contains no variants -> ref
                    calls2_c2, ref_c2, calls2_str_c2, ref_str_c2, 
                    calls2_ref_ptrs, ref_calls2_ptrs, calls2_vars, calls2_vars,
                    sc->calls2_beg_idx[sc_idx], 0,
                    sc->calls2_end_idx[sc_idx], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth1_t1 = "", ref_t1 = "", truth1_str_t1 = "", ref_str_t1 = ""; 
            std::vector<int> truth1_ref_ptrs, ref_truth1_ptrs;
            generate_ptrs_strs( // truth1_vars[0:0] contains no variants -> ref
                    truth1_t1, ref_t1, truth1_str_t1, ref_str_t1, 
                    truth1_ref_ptrs, ref_truth1_ptrs, truth1_vars, truth1_vars,
                    sc->truth1_beg_idx[sc_idx], 0,
                    sc->truth1_end_idx[sc_idx], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth2_t2 = "", ref_t2 = "", truth2_str_t2 = "", ref_str_t2 = ""; 
            std::vector<int> truth2_ref_ptrs, ref_truth2_ptrs;
            generate_ptrs_strs( // truth2_vars[0:0] contains no variants -> ref
                    truth2_t2, ref_t2, truth2_str_t2, ref_str_t2, 
                    truth2_ref_ptrs, ref_truth2_ptrs, truth2_vars, truth2_vars,
                    sc->truth2_beg_idx[sc_idx], 0,
                    sc->truth2_end_idx[sc_idx], 0,
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );

            // precision/recall alignment
            std::vector<int> pr_score(4), pr_calls_ref_end(4);
            std::vector< std::vector< std::vector<int> > > pr_offs(8), pr_ptrs;
            calc_prec_recall_aln(
                    calls1_c1, calls2_c2, truth1_t1, truth2_t2, ref_c1,
                    calls1_ref_ptrs, ref_calls1_ptrs, 
                    calls2_ref_ptrs, ref_calls2_ptrs,
                    pr_score, pr_ptrs, pr_calls_ref_end
            );

            // calculate optimal global phasing
            int phase = store_phase(clusterdata_ptr, ctg, pr_score);

            // EDIT DISTANCE
            // generate pointers and strings
            std::string calls1 = "", calls2 = "", calls1_str = "", calls2_str = ""; 
            std::vector<int> calls1_calls2_ptrs, calls2_calls1_ptrs;
            generate_ptrs_strs(
                    calls1, calls2, calls1_str, calls2_str, 
                    calls1_calls2_ptrs, calls2_calls1_ptrs, calls1_vars, calls2_vars,
                    sc->calls1_beg_idx[sc_idx], sc->calls2_beg_idx[sc_idx],
                    sc->calls1_end_idx[sc_idx], sc->calls2_end_idx[sc_idx],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth1 = "", truth2 = "", truth1_str = "", truth2_str = ""; 
            std::vector<int> truth1_truth2_ptrs, truth2_truth1_ptrs;
            generate_ptrs_strs(
                    truth1, truth2, truth1_str, truth2_str, 
                    truth1_truth2_ptrs, truth2_truth1_ptrs, truth1_vars, truth2_vars,
                    sc->truth1_beg_idx[sc_idx], sc->truth2_beg_idx[sc_idx],
                    sc->truth1_end_idx[sc_idx], sc->truth2_end_idx[sc_idx],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );

            // edit distance alignment
            std::vector<int> s(4);
            std::vector< std::vector< std::vector<int> > > offs(4), ptrs(4);
            calc_edit_dist_aln(calls1, calls2, truth1, truth2, s, offs, ptrs);

            // update total distance
            int orig_phase_dist = s[CALLS1_TRUTH1] + s[CALLS2_TRUTH2];
            int swap_phase_dist = s[CALLS2_TRUTH1] + s[CALLS1_TRUTH2];
            int dist = std::min(orig_phase_dist, swap_phase_dist);
            distance += dist;

            // calculate precision/recall from alignment
            color_ptrs(pr_ptrs, calls1_ref_ptrs, ref_calls1_ptrs, 
                    calls2_ref_ptrs, ref_calls2_ptrs, pr_calls_ref_end);
            calc_prec_recall(
                    clusterdata_ptr, sc_idx, ctg,
                    calls1, calls2, truth1, truth2, ref_c1,
                    calls1_ref_ptrs, ref_calls1_ptrs, 
                    calls2_ref_ptrs, ref_calls2_ptrs,
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    pr_ptrs, pr_calls_ref_end, phase, dist
            );

            // DEBUG PRINTING
            std::vector<std::string> calls {calls1, calls1, calls2, calls2};
            std::vector<std::string> truth {truth1, truth2, truth1, truth2};
            std::vector<int> calls_lens = 
                    {int(calls1.size()), int(calls1.size()), int(calls2.size()), int(calls2.size())};
            std::vector<int> truth_lens = 
                    {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};
            if (g.print_verbosity >= 1 && dist) {
                // print cluster info
                printf("\n\nCALLS1: %d clusters\n", sc->calls1_end_idx[sc_idx] - sc->calls1_beg_idx[sc_idx]);
                for(int i = sc->calls1_beg_idx[sc_idx]; i < sc->calls1_end_idx[sc_idx]; i++) {
                    printf("\tGroup %d: %d variants (%d-%d)\n", i, 
                            calls1_vars->clusters[i+1]-calls1_vars->clusters[i],
                            calls1_vars->clusters[i], calls1_vars->clusters[i+1]);
                    for(int j = calls1_vars->clusters[i]; j < calls1_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), calls1_vars->poss[j], 
                                calls1_vars->refs[j].size() ? calls1_vars->refs[j].data() : "_", 
                                calls1_vars->alts[j].size() ? calls1_vars->alts[j].data() : "_");
                    }
                }
                printf("CALLS2: %d clusters\n", sc->calls2_end_idx[sc_idx] - sc->calls2_beg_idx[sc_idx]);
                for(int i = sc->calls2_beg_idx[sc_idx]; i < sc->calls2_end_idx[sc_idx]; i++) {
                    printf("\tGroup %d: %d variants (%d-%d)\n", i, 
                            calls2_vars->clusters[i+1]-calls2_vars->clusters[i],
                            calls2_vars->clusters[i], calls2_vars->clusters[i+1]);
                    for(int j = calls2_vars->clusters[i]; j < calls2_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), calls2_vars->poss[j], 
                                calls2_vars->refs[j].size() ? calls2_vars->refs[j].data() : "_", 
                                calls2_vars->alts[j].size() ? calls2_vars->alts[j].data() : "_");
                    }
                }
                printf("TRUTH1: %d clusters\n", sc->truth1_end_idx[sc_idx] - sc->truth1_beg_idx[sc_idx]);
                for(int i = sc->truth1_beg_idx[sc_idx]; i < sc->truth1_end_idx[sc_idx]; i++) {
                    printf("\tGroup %d: %d variants (%d-%d)\n", i, 
                            truth1_vars->clusters[i+1]-truth1_vars->clusters[i],
                            truth1_vars->clusters[i], truth1_vars->clusters[i+1]);
                    for(int j = truth1_vars->clusters[i]; j < truth1_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), truth1_vars->poss[j], 
                                truth1_vars->refs[j].size() ? truth1_vars->refs[j].data() : "_", 
                                truth1_vars->alts[j].size() ? truth1_vars->alts[j].data() : "_");
                    }
                }
                printf("TRUTH2: %d clusters\n", sc->truth2_end_idx[sc_idx] - sc->truth2_beg_idx[sc_idx]);
                for(int i = sc->truth2_beg_idx[sc_idx]; i < sc->truth2_end_idx[sc_idx]; i++) {
                    printf("\tGroup %d: %d variants (%d-%d)\n", i, 
                            truth2_vars->clusters[i+1]-truth2_vars->clusters[i],
                            truth2_vars->clusters[i], truth2_vars->clusters[i+1]);
                    for(int j = truth2_vars->clusters[i]; j < truth2_vars->clusters[i+1]; j++) {
                        printf("\t\t%s %d\t%s\t%s\n", ctg.data(), truth2_vars->poss[j], 
                                truth2_vars->refs[j].size() ? truth2_vars->refs[j].data() : "_", 
                                truth2_vars->alts[j].size() ? truth2_vars->alts[j].data() : "_");
                    }
                }

                // generate ref string
                std::string ref_str = clusterdata_ptr->ref->fasta.at(ctg).
                    substr(sc->begs[sc_idx], sc->ends[sc_idx] - sc->begs[sc_idx]);

                // DEBUG STR/PTR PRINTING
                printf("ORIG:      %s\n", ref_str.data());
                printf("CALLS1:    %s\n", calls1_str.data());
                printf("CALLS2:    %s\n", calls2_str.data());
                printf("TRUTH1:    %s\n", truth1_str.data());
                printf("TRUTH2:    %s\n", truth2_str.data());
                printf("Edit Distance: %d\n", dist);

                printf("ORIG_C1:   %s\n", ref_str_c1.data());
                printf("CALLS1_C1: %s\n", calls1_str_c1.data());
                printf("CALLS1_REF: ");
                for(size_t i = 0; i < calls1_ref_ptrs.size(); i++) {
                    printf("%d ", calls1_ref_ptrs[i]);
                } printf("\n");
                printf("REF_CALLS1: ");
                for(size_t i = 0; i < ref_calls1_ptrs.size(); i++) {
                    printf("%d ", ref_calls1_ptrs[i]);
                } printf("\n");

                printf("ORIG_C2:   %s\n", ref_str_c2.data());
                printf("CALLS2_C2: %s\n", calls2_str_c2.data());
                printf("CALLS2_REF: ");
                for(size_t i = 0; i < calls2_ref_ptrs.size(); i++) {
                    printf("%d ", calls2_ref_ptrs[i]);
                } printf("\n");
                printf("REF_CALLS2: ");
                for(size_t i = 0; i < ref_calls2_ptrs.size(); i++) {
                    printf("%d ", ref_calls2_ptrs[i]);
                } printf("\n");

            }

            // MORE DEBUG PRINTING
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

        } // each cluster
    } // each contig
    INFO("Total edit distance: %d", distance);
}
