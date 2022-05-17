#include <string>
#include <vector>

#include "dist.h"
#include "print.h"


int edit_dist_realign(const vcfData* vcf, const fastaData* const ref) {

    // iterate over each haplotype
    int groups = 0;
    int new_ed_groups = 0;
    int old_ed = 0;
    int new_ed = 0;
    for (int h = 0; h < 2; h++) {

        // iterate over each contig
        for (auto itr = vcf->hapcalls[h].begin(); 
                itr != vcf->hapcalls[h].end(); itr++) {
            std::string ctg = itr->first;
            variantCalls vars = itr->second;
            if (vars.poss.size() == 0) continue;

            // iterate over each group of variants
            for (size_t var_grp = 0; var_grp < vars.gaps.size()-1; var_grp++) {
                int beg_idx = vars.gaps[var_grp];
                int end_idx = vars.gaps[var_grp+1];
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
                std::string ref_str = "";
                std::string alt_str = "";
                std::string alt = "";
                for (int ref_pos = beg; ref_pos < end;) {
                    if (ref_pos == vars.poss[var]) { // in variant
                        switch (vars.types[var]) {
                            case TYPE_INS:
                                alt += vars.alts[var];
                                alt_str += GREEN(vars.alts[var]);
                                ref_str += std::string(vars.alts[var].size(), ' ');
                                break;
                            case TYPE_DEL:
                                alt_str += std::string(vars.refs[var].size(), ' ');
                                ref_str += RED(vars.refs[var]);
                                ref_pos += vars.refs[var].size();
                                break;
                            case TYPE_SUB:
                                alt += vars.alts[var];
                                alt_str += " " + GREEN(vars.alts[var]);
                                ref_str += RED(vars.refs[var]) + " ";
                                ref_pos++;
                                break;
                            case TYPE_GRP:
                                alt += vars.alts[var];
                                alt_str += std::string(vars.refs[var].size(), ' ') 
                                    + GREEN(vars.alts[var]);
                                ref_str += RED(vars.refs[var]) + std::string(
                                        vars.alts[var].size(), ' ');
                                ref_pos += vars.refs[var].size();
                                break;
                        }
                        var++; // next variant
                    }
                    else { // match
                        try {
                            alt += ref->fasta.at(ctg)[ref_pos];
                            ref_str += ref->fasta.at(ctg)[ref_pos];
                            alt_str += ref->fasta.at(ctg)[ref_pos];
                            ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // do alignment
                int s = 0;
                int reflen = end-beg;
                int altlen = reflen + inss - dels;
                std::vector< std::vector<int> > diags, offsets, ptrs;
                diags.push_back(std::vector<int>(1,0));
                offsets.push_back(std::vector<int>(1,-1));
                ptrs.push_back(std::vector<int>(1,0));
                bool done = false;
                while (true) {

                    // EXTEND WAVEFRONT
                    for (int d = 0; d < s+1; d++) {
                        int max_offset = std::min(
                                reflen - diags[s][d], altlen) - 1;
                        int offset = offsets[s][d];
                        while (offset < max_offset && alt[offset+1] == 
                                ref->fasta.at(ctg)[diags[s][d]+offset+beg+1]) {
                            offset++;
                        }
                        offsets[s][d] = offset;
                        if (offset == altlen-1 && 
                                offset+diags[s][d] == reflen-1)
                        { done = true; break; }
                    }
                    if (done) break;


                    // NEXT WAVEFRONT
                    
                    // add wavefront, fill edge cells
                    diags.push_back(std::vector<int>(s+2));
                    offsets.push_back(std::vector<int>(s+2));
                    ptrs.push_back(std::vector<int>(s+2));
                    diags[s+1][0] = diags[s][0] - 1;
                    offsets[s+1][0] = offsets[s][0] + 1;
                    ptrs[s+1][0] = PTR_UP;
                    diags[s+1][s+1] = diags[s][s] + 1;
                    offsets[s+1][s+1] = offsets[s][s];
                    ptrs[s+1][s+1] = PTR_LEFT;

                    // central cells
                    for (int d = 1; d <= s; d++) {
                        if (offsets[s][d-1] >= offsets[s][d]+1) {
                            diags[s+1][d] = diags[s][d-1] + 1;
                            offsets[s+1][d] = offsets[s][d-1];
                            ptrs[s+1][d] = PTR_LEFT;
                        } else {
                            diags[s+1][d] = diags[s][d] - 1;
                            offsets[s+1][d] = offsets[s][d]+1;
                            ptrs[s+1][d] = PTR_UP;
                        }
                    }
                    ++s;
                }

                // DEBUG PRINT

                if (g.print_verbosity >= 2) {
                    printf("\noffsets:\n");
                    printf("s= ");
                    for(int i = 0; i <= s; i++) printf("%2i ", i);
                    printf("\n ");
                    for(int i = 0; i <= s; i++) printf("  /");
                    printf("\n");
                    for(int r = 0; r <= s; r++) {
                        for(int c = 0; c <= s-r; c++) {
                            printf(" %2d", offsets[r+c][c]);
                        }
                        printf("\n");
                    }

                    // create array
                    std::vector< std::vector<char> > ptr_str;
                    for (int i = 0; i < altlen; i++)
                        ptr_str.push_back(std::vector<char>(reflen, '.'));

                    // modify array with pointers
                    int altpos, refpos;
                    for (int si = 0; si <= s; si++) {
                        for(int di = 0; di <= si; di++) {
                            if (di == 0) {
                                if (si == 0) {
                                    altpos = 0;
                                    refpos = 0;
                                } else {
                                    altpos = offsets[si-1][di] + 1;
                                    refpos = diags[si-1][di] + offsets[si-1][di];
                                    if (altpos < altlen && refpos < reflen)
                                        ptr_str[altpos][refpos] = '|';
                                }
                            } 
                            else if (di == si) {
                                altpos = offsets[si-1][di-1];
                                refpos = diags[si-1][di-1] + offsets[si-1][di-1] + 1;
                                if (altpos < altlen && refpos < reflen)
                                    ptr_str[altpos][refpos] = '-';
                            } 
                            else if (offsets[si-1][di-1] > offsets[si-1][di]+1) {
                                altpos = offsets[si-1][di-1];
                                refpos = diags[si-1][di-1] + offsets[si-1][di-1] + 1;
                                if (altpos < altlen && refpos < reflen)
                                    ptr_str[altpos][refpos] = '-';
                            } 
                            else {
                                altpos = offsets[si-1][di] + 1;
                                refpos = diags[si-1][di] + offsets[si-1][di];
                                if (altpos < altlen && refpos < reflen)
                                    ptr_str[altpos][refpos] = '|';
                            }
                            while (altpos < altlen-1 && refpos < reflen-1 && 
                                    alt[++altpos] == ref->fasta.at(ctg)[beg + ++refpos]) {
                                ptr_str[altpos][refpos] = '\\';
                            }
                        }
                    }
                    // NOTE: first base will ALWAYS match, due to grabbing 
                    // previous ref base before variant. This is required 
                    // for correctness of algorithm
                    ptr_str[0][0] = '\\';

                    // print array
                    for (int i = -1; i < altlen; i++) {
                        for (int j = -1; j < reflen; j++) {
                            if (i < 0 && j < 0) {
                                printf("\n  ");
                            }
                            else if (i < 0) {
                                printf("%c", ref->fasta.at(ctg)[beg+j]);
                            } else if (j < 0) {
                                printf("\n%c ", alt[i]);
                            } else {
                                printf("%c", ptr_str[i][j]);
                            }
                        }
                    }
                }

                // BACKTRACK
                
                // init
                std::vector<int> cig(reflen+altlen);
                int cig_ptr = reflen + altlen - 1;
                int score = s;
                int diag = (reflen-altlen+score) / 2; // idx of last cell in WF
                int ptr;

                while (score > 0) {

                    // find previous wavefront
                    int off = offsets[score][diag];
                    int prev_off, up_off, left_off;
                    if (diag == 0) { // left edge, must go up
                        if (score == 0) {
                            ERROR("score should not be zero.");
                        }
                        up_off = offsets[score-1][diag]+1;
                        prev_off = up_off;
                        ptr = PTR_UP;

                    } else if (diag == score) { // right edge, must go left
                        left_off = offsets[score-1][diag-1];
                        prev_off = left_off;
                        ptr = PTR_LEFT;

                    } else { // get predecessor
                        left_off = offsets[score-1][diag-1];
                        up_off = offsets[score-1][diag]+1;
                        if (left_off > up_off) {
                            prev_off = left_off;
                            ptr = PTR_LEFT;
                        } else {
                            prev_off = up_off;
                            ptr = PTR_UP;
                        }
                    }

                    // slide up diagonally
                    while (off > prev_off) {
                        cig[cig_ptr--] = PTR_DIAG;
                        cig[cig_ptr--] = PTR_DIAG;
                        off--;
                    }

                    // go to previous wavefront
                    cig[cig_ptr--] = ptr;
                    score--;
                    switch(ptr) {
                    case PTR_LEFT:
                        diag--;
                        break;
                    case PTR_UP:
                        off--;
                        break;
                    default:
                        ERROR("Pointer type unexpected: ptr=%i", ptr);
                        break;
                    }

                }
                // slide up diagonally to end
                for (int i = 0; i <= offsets[0][0]; i++) {
                    cig[cig_ptr--] = PTR_DIAG;
                    cig[cig_ptr--] = PTR_DIAG;
                }

                // get ref/alt strings for printing
                int ref_ptr = 0;
                int alt_ptr = 0;
                int new_inss = 0;
                int new_dels = 0;
                std::string new_ref_str, new_alt_str;
                for(size_t i = 0; i < cig.size(); i++) {
                    switch(cig[i]) {
                        case PTR_DIAG:
                            new_ref_str += ref->fasta.at(ctg)[beg+ref_ptr++];
                            new_alt_str += alt[alt_ptr++];
                            i++;
                            break;
                        case PTR_UP:
                            new_ref_str += " ";
                            new_alt_str += GREEN(alt[alt_ptr++]);
                            new_inss++;
                            break;
                        case PTR_LEFT:
                            new_ref_str += RED(ref->fasta.at(ctg)[beg+ref_ptr++]);
                            new_alt_str += " ";
                            new_dels++;
                            break;
                    }
                }

                if (g.print_verbosity >= 2) {
                    printf("\nPTRS: ");
                    for(size_t i = 0; i < cig.size(); i++) {
                        printf("%i ", cig[i]);
                    }

                    printf("\nCIGAR: ");
                    for(size_t i = 0; i < cig.size(); i++) {
                        switch(cig[i]) {
                            case PTR_DIAG:
                                if (cig[++i] != PTR_DIAG)
                                    ERROR("cig should be PTR_DIAG");
                                printf("M"); break;
                            case PTR_UP:
                                printf("I"); break;
                            case PTR_LEFT:
                                printf("D"); break;
                        }
                    }
                    printf("\n");
                }

                // print group info, and old/new alignments
                if (g.print_verbosity >= 1) {
                    if (subs*2 + inss + dels != s) {
                        const char* var_ref;
                        const char* var_alt;
                        printf("\n\n  Group %i: %d variants, %s:%d-%d\n\n"
                                "    old edit distance: %d (%dS %dI %dD)\n",
                                int(var_grp), end_idx-beg_idx, ctg.data(), 
                                beg, end, subs*2 + inss + dels, subs, inss, dels);
                        for (int var = beg_idx; var < end_idx; var++) {
                            if (vars.refs[var].size() == 0) var_ref = "-"; 
                            else var_ref = vars.refs[var].data();
                            if (vars.alts[var].size() == 0) var_alt = "-"; 
                            else var_alt = vars.alts[var].data();
                            printf("      %s:%i hap%i %s\t%s\t%s\n", 
                                    ctg.data(), vars.poss[var], h, 
                                    type_strs[vars.types[var]].data(), 
                                    var_ref, var_alt);
                        }
                        printf("      REF: %s\n      ALT: %s\n", 
                                ref_str.data(), alt_str.data());
                        printf("\n    new edit distance: %i (%iI %iD)\n", 
                                s, new_inss, new_dels);
                        printf("      REF: %s\n      ALT: %s\n", 
                                new_ref_str.data(), new_alt_str.data());
                    }
                }

                // update counters
                old_ed += subs*2 + inss + dels;
                new_ed += s;
                if (subs*2 + inss + dels > s) new_ed_groups++;
                groups++;

            } // group
        } // contig
    } // hap
    INFO("Edit dist reduced in %i of %i groups, from %i to %i.", 
            new_ed_groups, groups, old_ed, new_ed);
    return 0;

}
