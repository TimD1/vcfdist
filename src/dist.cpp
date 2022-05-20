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

/******************************************************************************/

int edit_dist(vcfData* calls, vcfData* truth, fastaData* ref) {

    // iterate over each haplotype
    int distance = 0;
    for (int call_hap = 0; call_hap < 2; call_hap++) {

        // iterate over each contig
        for (auto itr = calls->hapcalls[call_hap].begin(); 
                itr != calls->hapcalls[call_hap].end(); itr++) {

            // initialize variant call info
            std::string ctg = itr->first;
            variantCalls* call_vars = &itr->second;
            if (call_vars->poss.size() == 0) continue;
            variantCalls* hap1_vars;
            variantCalls* hap2_vars;
            try {
                hap1_vars = &truth->hapcalls[0][ctg];
                hap2_vars = &truth->hapcalls[1][ctg];
            } catch (const std::exception & e) {
                ERROR("truth VCF does not contain contig '%s'", ctg.data());
                return EXIT_FAILURE;
            }
            size_t call_clust_beg_idx = 0;
            size_t hap1_clust_beg_idx = 0;
            size_t hap2_clust_beg_idx = 0;

            // for each cluster of variants (merge calls and truth haps)
            while (call_clust_beg_idx < call_vars->gaps.size()-1 ||
                   hap1_clust_beg_idx < hap1_vars->gaps.size()-1 ||
                   hap2_clust_beg_idx < hap2_vars->gaps.size()-1) {

                // set all start positions
                size_t call_clust_end_idx = call_clust_beg_idx;
                size_t hap1_clust_end_idx = hap1_clust_beg_idx;
                size_t hap2_clust_end_idx = hap2_clust_beg_idx;
                int next_call_beg_pos, next_hap1_beg_pos, next_hap2_beg_pos;
                if (call_clust_end_idx >= call_vars->gaps.size()-1)
                    next_call_beg_pos = std::numeric_limits<int>::max();
                else
                    next_call_beg_pos = call_vars->poss[call_vars->gaps[call_clust_end_idx]]-1;
                if (hap1_clust_end_idx >= hap1_vars->gaps.size()-1)
                    next_hap1_beg_pos = std::numeric_limits<int>::max();
                else
                    next_hap1_beg_pos = hap1_vars->poss[hap1_vars->gaps[hap1_clust_end_idx]]-1;
                if (hap2_clust_end_idx >= hap2_vars->gaps.size()-1)
                    next_hap2_beg_pos = std::numeric_limits<int>::max();
                else
                    next_hap2_beg_pos = hap2_vars->poss[hap2_vars->gaps[hap2_clust_end_idx]]-1;

                // initialize cluster merging with first to start
                int curr_end_pos = 0;
                if (next_call_beg_pos < next_hap1_beg_pos && 
                        next_call_beg_pos < next_hap2_beg_pos) { // call var first
                    call_clust_end_idx += 1;
                    curr_end_pos = call_vars->poss[call_vars->gaps[call_clust_end_idx]-1] + 
                            call_vars->rlens[call_vars->gaps[call_clust_end_idx]-1] + 1;
                    if (call_clust_end_idx >= call_vars->gaps.size()-1)
                        next_call_beg_pos = std::numeric_limits<int>::max();
                    else
                        next_call_beg_pos = call_vars->poss[call_vars->gaps[call_clust_end_idx]]-1;
                } 
                else if (next_hap1_beg_pos < next_call_beg_pos && 
                        next_hap1_beg_pos < next_hap2_beg_pos) { // hap1 var first
                    hap1_clust_end_idx += 1;
                    curr_end_pos = hap1_vars->poss[hap1_vars->gaps[hap1_clust_end_idx]-1] + 
                            hap1_vars->rlens[hap1_vars->gaps[hap1_clust_end_idx]-1] + 1;
                    if (hap1_clust_end_idx >= hap1_vars->gaps.size()-1)
                        next_hap1_beg_pos = std::numeric_limits<int>::max();
                    else
                        next_hap1_beg_pos = hap1_vars->poss[hap1_vars->gaps[hap1_clust_end_idx]]-1;
                } 
                else { // hap2 var first
                    hap2_clust_end_idx += 1;
                    curr_end_pos = hap2_vars->poss[hap2_vars->gaps[hap2_clust_end_idx]-1] + 
                            hap2_vars->rlens[hap2_vars->gaps[hap2_clust_end_idx]-1] + 1;
                    if (hap2_clust_end_idx >= hap2_vars->gaps.size()-1)
                        next_hap2_beg_pos = std::numeric_limits<int>::max();
                    else
                        next_hap2_beg_pos = hap2_vars->poss[hap2_vars->gaps[hap2_clust_end_idx]]-1;
                }

                // keep expanding cluster while possible
                bool just_merged = true;
                bool any_merged = false;
                while (just_merged) {
                    just_merged = false;
                    while (next_hap1_beg_pos < curr_end_pos + g.gap) {
                        hap1_clust_end_idx += 1;
                        curr_end_pos = hap1_vars->poss[hap1_vars->gaps[hap1_clust_end_idx]-1] + 
                                hap1_vars->rlens[hap1_vars->gaps[hap1_clust_end_idx]-1] + 1;
                        if (hap1_clust_end_idx >= hap1_vars->gaps.size()-1)
                            next_hap1_beg_pos = std::numeric_limits<int>::max();
                        else
                            next_hap1_beg_pos = hap1_vars->poss[hap1_vars->gaps[hap1_clust_end_idx]]-1;
                        just_merged = true;
                    }
                    while (next_hap2_beg_pos < curr_end_pos + g.gap) {
                        hap2_clust_end_idx += 1;
                        curr_end_pos = hap2_vars->poss[hap2_vars->gaps[hap2_clust_end_idx]-1] + 
                                hap2_vars->rlens[hap2_vars->gaps[hap2_clust_end_idx]-1] + 1;
                        if (hap2_clust_end_idx >= hap2_vars->gaps.size()-1)
                            next_hap2_beg_pos = std::numeric_limits<int>::max();
                        else
                            next_hap2_beg_pos = hap2_vars->poss[hap2_vars->gaps[hap2_clust_end_idx]]-1;
                        just_merged = true;
                    }
                    while (next_call_beg_pos < curr_end_pos + g.gap) {
                        call_clust_end_idx += 1;
                        curr_end_pos = call_vars->poss[call_vars->gaps[call_clust_end_idx]-1] + 
                                call_vars->rlens[call_vars->gaps[call_clust_end_idx]-1] + 1;
                        if (call_clust_end_idx >= call_vars->gaps.size()-1)
                            next_call_beg_pos = std::numeric_limits<int>::max();
                        else
                            next_call_beg_pos = call_vars->poss[call_vars->gaps[call_clust_end_idx]]-1;
                        just_merged = true;
                    }
                    if (just_merged) any_merged = true;
                }

                // get supercluster start/end positions (allowing empty haps)
                int beg_pos = std::numeric_limits<int>::max();
                int end_pos = -1;
                if (call_clust_end_idx - call_clust_beg_idx) { // call vars present
                    beg_pos = std::min(beg_pos, 
                            call_vars->poss[call_vars->gaps[call_clust_beg_idx]]-1);
                    end_pos = std::max(end_pos, 
                            call_vars->poss[call_vars->gaps[call_clust_end_idx]-1] + 
                            call_vars->rlens[call_vars->gaps[call_clust_end_idx]-1]+1);
                }
                if (hap1_clust_end_idx - hap1_clust_beg_idx) { // hap1 vars present
                    beg_pos = std::min(beg_pos, 
                            hap1_vars->poss[hap1_vars->gaps[hap1_clust_beg_idx]]-1);
                    end_pos = std::max(end_pos, 
                            hap1_vars->poss[hap1_vars->gaps[hap1_clust_end_idx]-1] + 
                            hap1_vars->rlens[hap1_vars->gaps[hap1_clust_end_idx]-1]+1);
                }
                if (hap2_clust_end_idx - hap2_clust_beg_idx) { // hap2 vars present
                    beg_pos = std::min(beg_pos, 
                            hap2_vars->poss[hap2_vars->gaps[hap2_clust_beg_idx]]-1);
                    end_pos = std::max(end_pos, 
                            hap2_vars->poss[hap2_vars->gaps[hap2_clust_end_idx]-1] + 
                            hap2_vars->rlens[hap2_vars->gaps[hap2_clust_end_idx]-1]+1);
                }

                // generate call string
                int call_var_idx = call_vars->gaps[call_clust_beg_idx];
                std::string call = "", call_str = "", ref_str = "";
                for (int ref_pos = beg_pos; ref_pos < end_pos; ) {
                    if (ref_pos == call_vars->poss[call_var_idx]) { // in call variant
                        switch (call_vars->types[call_var_idx]) {
                            case TYPE_INS:
                                call += call_vars->alts[call_var_idx];
                                call_str += GREEN(call_vars->alts[call_var_idx]);
                                ref_str += std::string(call_vars->alts[call_var_idx].size(), ' ');
                                break;
                            case TYPE_DEL:
                                call_str += std::string(call_vars->refs[call_var_idx].size(), ' ');
                                ref_pos += call_vars->refs[call_var_idx].size();
                                ref_str += RED(call_vars->refs[call_var_idx]);
                                break;
                            case TYPE_SUB:
                                call += call_vars->alts[call_var_idx];
                                call_str += GREEN(call_vars->alts[call_var_idx]);
                                ref_str += RED(call_vars->refs[call_var_idx]);
                                ref_pos++;
                                break;
                            case TYPE_GRP:
                                call += call_vars->alts[call_var_idx];
                                call_str += std::string(call_vars->refs[call_var_idx].size(), ' ') 
                                    + GREEN(call_vars->alts[call_var_idx]);
                                ref_str += RED(call_vars->refs[call_var_idx]) + std::string(
                                        call_vars->alts[call_var_idx].size(), ' ');
                                ref_pos += call_vars->refs[call_var_idx].size();
                                break;
                        }
                        call_var_idx++; // next variant
                    }
                    else { // match
                        try {
                            call += ref->fasta.at(ctg)[ref_pos];
                            call_str += ref->fasta.at(ctg)[ref_pos];
                            ref_str += ref->fasta.at(ctg)[ref_pos];
                            ref_pos++;
                        } catch (const std::out_of_range & e) {
                            ERROR("contig %s not present in reference FASTA",
                                    ctg.data());
                            exit(1);
                        }
                    }
                }

                // generate hap1 and hap2 strings and pointers
                int hap1_var_idx = hap1_vars->gaps[hap1_clust_beg_idx];
                int hap2_var_idx = hap2_vars->gaps[hap2_clust_beg_idx];
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

                // do alignment
                int s = 0;
                int call_len = call.size();
                int hap1_len = hap1.size();
                int hap2_len = hap2.size();
                int hap1_mat_len = call.size() + hap1.size() - 1;
                int hap2_mat_len = call.size() + hap2.size() - 1;
                std::vector< std::vector<int> > offsets1, ptrs1;
                offsets1.push_back(std::vector<int>(hap1_mat_len,-2)); // init all invalid
                /* ptrs1.push_back(std::vector<int>(hap1_mat_len,PTR_NONE)); */
                std::vector< std::vector<int> > offsets2, ptrs2;
                offsets2.push_back(std::vector<int>(hap2_mat_len,-2)); // init all invalid
                /* ptrs2.push_back(std::vector<int>(hap2_mat_len,PTR_NONE)); */
                offsets1[0][call_len-1] = -1; // only start from top left diag
                offsets2[0][call_len-1] = -1;
                bool done = false;
                bool extend = true;


                while (true) {

                    // EXTEND WAVEFRONT
                    extend = true;
                    while (extend) {

                        // extend hap1, allow hap1 -> hap2 transition
                        for (int d1 = 0; d1 < hap1_mat_len; d1++) {
                            int off1 = offsets1[s][d1];
                            int diag1 = d1 + 1 - call_len;

                            // don't allow starting from untouched cells
                            if (off1 == -2) continue;

                            // check that it's within matrix
                            if (diag1 + off1 + 1 < 0) continue;
                            if (off1 > call_len - 1) continue;
                            if (diag1 + off1 > hap1_len - 1) continue;

                            // extend
                            while (off1 < call_len - 1 && 
                                   diag1 + off1 < hap1_len - 1) {

                                // extend to other hap if possible
                                if (off1 + diag1 < hap1_len &&
                                    off1 + diag1 >= 0 &&
                                    hap1_ptrs[off1+diag1] >= 0) {
                                    int d2 = hap1_ptrs[off1+diag1] - off1 + call_len-1;
                                    if (off1 > offsets2[s][d2]) {
                                        /* printf("switch 1->2: off1=%d offsets2=%d s=%d d2=%d diag1=%d\n", */
                                        /*         off1, offsets2[s][d2], s, d2, diag1); */
                                        offsets2[s][d2] = off1;
                                    }
                                }

                                // extend on this hap if possible
                                if (call[off1+1] == hap1[diag1+off1+1])
                                    off1++;
                                else
                                    break;
                            }
                            offsets1[s][d1] = off1;

                            // finish if done
                            if (off1 == call_len - 1 && 
                                off1 + diag1 == hap1_len - 1)
                            { done = true; break; }

                        }
                        if (done) break;

                        // extend hap2, allow hap2->hap1 transition, redo if taken
                        extend = false;
                        for (int d2 = 0; d2 < hap2_mat_len; d2++) {
                            int off2 = offsets2[s][d2];
                            int diag2 = d2 + 1 - call_len;

                            // don't allow starting from untouched cells
                            if (off2 == -2) continue;

                            // check that it's within matrix
                            if (diag2 + off2 + 1 < 0) continue;
                            if (off2 > call_len - 1) continue;
                            if (diag2 + off2 > hap2_len - 1) continue;

                            // extend
                            while (off2 < call_len - 1 && 
                                   diag2 + off2 < hap2_len - 1) {

                                // extend to other hap if possible
                                if (off2 + diag2 < hap2_len &&
                                    off2 + diag2 >= 0 &&
                                    hap2_ptrs[off2+diag2] >= 0) {
                                    int d1 = hap2_ptrs[off2+diag2] - off2 + call_len-1;
                                    if (off2 > offsets1[s][d1]) {
                                        extend = true;
                                        /* printf("switch 2->1: off2=%d offsets1=%d s=%d d1=%d diag2=%d\n", */
                                        /*         off2, offsets1[s][d1], s, d1, diag2); */
                                        offsets1[s][d1] = off2;
                                    }
                                }

                                // extend on this hap if possible
                                if (call[off2+1] == hap2[diag2+off2+1])
                                    off2++;
                                else
                                    break;
                            }
                            offsets2[s][d2] = off2;

                            // finish if done
                            if (off2 == call_len - 1 && 
                                off2 + diag2 == hap2_len - 1)
                            { done = true; break; }

                        }
                        if (done) break;
                    }
                    if (done) break;


                    // NEXT WAVEFRONT
                    
                    // add wavefront, fill edge cells
                    offsets1.push_back(std::vector<int>(hap1_mat_len, -2));
                    offsets2.push_back(std::vector<int>(hap2_mat_len, -2));
                    // bottom left cells
                    if (s+1 == call_len-1) {
                        offsets1[s+1][0] = s+1;
                        offsets2[s+1][0] = s+1;
                    }
                    // top right cells
                    if (s+1 == hap1_mat_len-1)
                        offsets1[s+1][hap1_mat_len-1] = s+1;
                    if (s+1 == hap2_mat_len-1)
                        offsets2[s+1][hap2_mat_len-1] = s+1;

                    /* ptrs1.push_back(std::vector<int>(hap1_mat_len)); */
                    /* ptrs1[s+1][0] = PTR_UP; */
                    /* ptrs1[s+1][s+1] = PTR_LEFT; */

                    // central cells
                    for (int d1 = 1; d1 < hap1_mat_len-1; d1++) {
                        int offleft1 = offsets1[s][d1-1];
                        int offtop1 = (offsets1[s][d1+1] == -2) ? 
                            -2 : offsets1[s][d1+1]+1;
                        if (offleft1 >= offtop1) {
                            offsets1[s+1][d1] = offleft1;
                            /* ptrs1[s+1][d1] = PTR_LEFT; */
                        } else {
                            offsets1[s+1][d1] = offtop1;
                            /* ptrs1[s+1][d1] = PTR_UP; */
                        }
                    }
                    for (int d2 = 1; d2 < hap2_mat_len-1; d2++) {
                        int offleft2 = offsets2[s][d2-1];
                        int offtop2 = (offsets2[s][d2+1] == -2) ? 
                            -2 : offsets2[s][d2+1]+1;
                        if (offleft2 >= offtop2) {
                            offsets2[s+1][d2] = offleft2;
                            /* ptrs1[s+1][d2] = PTR_LEFT; */
                        } else {
                            offsets2[s+1][d2] = offtop2;
                            /* ptrs1[s+1][d2] = PTR_UP; */
                        }
                    }
                    ++s;
                }


                // PRINT RESULTS
                if (s > 0 && any_merged) 
                {

                    // print cluster info
                    printf("\n\nCALL: %zu groups\n", call_clust_end_idx-call_clust_beg_idx);
                    for(size_t i = call_clust_beg_idx; i < call_clust_end_idx; i++) {
                        printf("\tGroup %zu: %d variants\n", i, call_vars->gaps[i+1]-call_vars->gaps[i]);
                        for(int j = call_vars->gaps[i]; j < call_vars->gaps[i+1]; j++) {
                            printf("\t\t%s %d\t%s\t%s\n", ctg.data(), call_vars->poss[j], 
                                    call_vars->refs[j].size() ? call_vars->refs[j].data() : "_", 
                                    call_vars->alts[j].size() ? call_vars->alts[j].data() : "_");
                        }
                    }
                    printf("HAP1: %zu groups\n", hap1_clust_end_idx-hap1_clust_beg_idx);
                    for(size_t i = hap1_clust_beg_idx; i < hap1_clust_end_idx; i++) {
                        printf("\tGroup %zu: %d variants\n", i, hap1_vars->gaps[i+1]-hap1_vars->gaps[i]);
                        for(int j = hap1_vars->gaps[i]; j < hap1_vars->gaps[i+1]; j++) {
                            printf("\t\t%s %d\t%s\t%s\n", ctg.data(), hap1_vars->poss[j], 
                                    hap1_vars->refs[j].size() ? hap1_vars->refs[j].data() : "_", 
                                    hap1_vars->alts[j].size() ? hap1_vars->alts[j].data() : "_");
                        }
                    }
                    printf("HAP2: %zu groups\n", hap2_clust_end_idx-hap2_clust_beg_idx);
                    for(size_t i = hap2_clust_beg_idx; i < hap2_clust_end_idx; i++) {
                        printf("\tGroup %zu: %d variants\n", i, hap2_vars->gaps[i+1]-hap2_vars->gaps[i]);
                        for(int j = hap2_vars->gaps[i]; j < hap2_vars->gaps[i+1]; j++) {
                            printf("\t\t%s %d\t%s\t%s\n", ctg.data(), hap2_vars->poss[j], 
                                    hap2_vars->refs[j].size() ? hap2_vars->refs[j].data() : "_", 
                                    hap2_vars->alts[j].size() ? hap2_vars->alts[j].data() : "_");
                        }
                    }
                    
                    printf("ORIG: %s\n", ref_str.data());
                    printf("CALL: %s\n", call_str.data());
                    printf("HAP1: %s\n", hap1_str.data());
                    printf("HAP2: %s\n", hap2_str.data());
                    printf("distance: %d\n", s);

                    /* // DEBUG PRINT */
                    /* printf("\noffsets hap 1:\n"); */
                    /* printf("s= "); */
                    /* for(int i = 0; i < hap1_mat_len; i++) printf("%2i ", i); */
                    /* printf("\n "); */
                    /* for(int i = 0; i < hap1_mat_len; i++) printf("  /"); */
                    /* printf("\n"); */
                    /* for(int r = 0; r < s; r++) { */
                    /*     for(int c = 0; c < hap1_mat_len; c++) { */
                    /*         printf(" %2d", offsets1[r][c]); */
                    /*     } */
                    /*     printf("\n"); */
                    /* } */
                    /* printf("\noffsets hap 2:\n"); */
                    /* printf("s= "); */
                    /* for(int i = 0; i < hap2_mat_len; i++) printf("%2i ", i); */
                    /* printf("\n "); */
                    /* for(int i = 0; i < hap2_mat_len; i++) printf("  /"); */
                    /* printf("\n"); */
                    /* for(int r = 0; r < s; r++) { */
                    /*     for(int c = 0; c < hap2_mat_len; c++) { */
                    /*         printf(" %2d", offsets2[r][c]); */
                    /*     } */
                    /*     printf("\n"); */
                    /* } */

                    /* // create array */
                    /* std::vector< std::vector<char> > ptr_str; */
                    /* for (int i = 0; i < call_len; i++) */
                    /*     ptr_str.push_back(std::vector<char>(hap1_len, '.')); */

                    /* // modify array with pointers */
                    /* int call_pos, hap1_pos; */
                    /* for (int si = 0; si <= s; si++) { */
                    /*     for(int di = 0; di <= si; di++) { */
                    /*         if (di == 0) { */
                    /*             if (si == 0) { */
                    /*                 call_pos = 0; */
                    /*                 hap1_pos = 0; */
                    /*             } else { */
                    /*                 call_pos = offsets1[si-1][di] + 1; */
                    /*                 hap1_pos = diags1[si-1][di] + offsets1[si-1][di]; */
                    /*                 if (call_pos < call_len && hap1_pos < hap1_len) */
                    /*                     ptr_str[call_pos][hap1_pos] = '|'; */
                    /*             } */
                    /*         } */ 
                    /*         else if (di == si) { */
                    /*             call_pos = offsets1[si-1][di-1]; */
                    /*             hap1_pos = diags1[si-1][di-1] + offsets1[si-1][di-1] + 1; */
                    /*             if (call_pos < call_len && hap1_pos < hap1_len) */
                    /*                 ptr_str[call_pos][hap1_pos] = '-'; */
                    /*         } */ 
                    /*         else if (offsets1[si-1][di-1] > offsets1[si-1][di]+1) { */
                    /*             call_pos = offsets1[si-1][di-1]; */
                    /*             hap1_pos = diags1[si-1][di-1] + offsets1[si-1][di-1] + 1; */
                    /*             if (call_pos < call_len && hap1_pos < hap1_len) */
                    /*                 ptr_str[call_pos][hap1_pos] = '-'; */
                    /*         } */ 
                    /*         else { */
                    /*             call_pos = offsets1[si-1][di] + 1; */
                    /*             hap1_pos = diags1[si-1][di] + offsets1[si-1][di]; */
                    /*             if (call_pos < call_len && hap1_pos < hap1_len) */
                    /*                 ptr_str[call_pos][hap1_pos] = '|'; */
                    /*         } */
                    /*         while (call_pos < call_len-1 && hap1_pos < hap1_len-1 && */ 
                    /*                 call[++call_pos] == hap1[++hap1_pos]) { */
                    /*             ptr_str[call_pos][hap1_pos] = '\\'; */
                    /*         } */
                    /*     } */
                    /* } */
                    /* // NOTE: first base will ALWAYS match, due to grabbing */ 
                    /* // previous ref base before variant. This is required */ 
                    /* // for correctness of algorithm */
                    /* ptr_str[0][0] = '\\'; */

                    /* // print array */
                    /* for (int i = -1; i < call_len; i++) { */
                    /*     for (int j = -1; j < hap1_len; j++) { */
                    /*         if (i < 0 && j < 0) { */
                    /*             printf("\n  "); */
                    /*         } */
                    /*         else if (i < 0) { */
                    /*             printf("%c", hap1[j]); */
                    /*         } else if (j < 0) { */
                    /*             printf("\n%c ", call[i]); */
                    /*         } else { */
                    /*             printf("%c", ptr_str[i][j]); */
                    /*         } */
                    /*     } */
                    /* } */

                }

                // reset for next merged cluster
                call_clust_beg_idx = call_clust_end_idx;
                hap1_clust_beg_idx = hap1_clust_end_idx;
                hap2_clust_beg_idx = hap2_clust_end_idx;
                distance += s;
            }
            INFO("Total edit distance: %d", distance);
        }
    }

    return distance;
}



            /*     // BACKTRACK */
                
            /*     // init */
            /*     std::vector<int> cig(reflen+altlen); */
            /*     int cig_ptr = reflen + altlen - 1; */
            /*     int score = s; */
            /*     int diag = (reflen-altlen+score) / 2; // idx of last cell in WF */
            /*     int ptr; */

            /*     while (score > 0) { */

            /*         // find previous wavefront */
            /*         int off = offsets[score][diag]; */
            /*         int prev_off, up_off, left_off; */
            /*         if (diag == 0) { // left edge, must go up */
            /*             if (score == 0) { */
            /*                 ERROR("score should not be zero."); */
            /*             } */
            /*             up_off = offsets[score-1][diag]+1; */
            /*             prev_off = up_off; */
            /*             ptr = PTR_UP; */

            /*         } else if (diag == score) { // right edge, must go left */
            /*             left_off = offsets[score-1][diag-1]; */
            /*             prev_off = left_off; */
            /*             ptr = PTR_LEFT; */

            /*         } else { // get predecessor */
            /*             left_off = offsets[score-1][diag-1]; */
            /*             up_off = offsets[score-1][diag]+1; */
            /*             if (left_off > up_off) { */
            /*                 prev_off = left_off; */
            /*                 ptr = PTR_LEFT; */
            /*             } else { */
            /*                 prev_off = up_off; */
            /*                 ptr = PTR_UP; */
            /*             } */
            /*         } */

            /*         // slide up diagonally */
            /*         while (off > prev_off) { */
            /*             cig[cig_ptr--] = PTR_DIAG; */
            /*             cig[cig_ptr--] = PTR_DIAG; */
            /*             off--; */
            /*         } */

            /*         // go to previous wavefront */
            /*         cig[cig_ptr--] = ptr; */
            /*         score--; */
            /*         switch(ptr) { */
            /*         case PTR_LEFT: */
            /*             diag--; */
            /*             break; */
            /*         case PTR_UP: */
            /*             off--; */
            /*             break; */
            /*         default: */
            /*             ERROR("Pointer type unexpected: ptr=%i", ptr); */
            /*             break; */
            /*         } */

            /*     } */
            /*     // slide up diagonally to end */
            /*     for (int i = 0; i <= offsets[0][0]; i++) { */
            /*         cig[cig_ptr--] = PTR_DIAG; */
            /*         cig[cig_ptr--] = PTR_DIAG; */
            /*     } */

            /*     // get ref/alt strings for printing */
            /*     int ref_ptr = 0; */
            /*     int alt_ptr = 0; */
            /*     int new_inss = 0; */
            /*     int new_dels = 0; */
            /*     std::string new_ref_str, new_alt_str; */
            /*     for(size_t i = 0; i < cig.size(); i++) { */
            /*         switch(cig[i]) { */
            /*             case PTR_DIAG: */
            /*                 new_ref_str += ref->fasta.at(ctg)[beg+ref_ptr++]; */
            /*                 new_alt_str += alt[alt_ptr++]; */
            /*                 i++; */
            /*                 break; */
            /*             case PTR_UP: */
            /*                 new_ref_str += " "; */
            /*                 new_alt_str += GREEN(alt[alt_ptr++]); */
            /*                 new_inss++; */
            /*                 break; */
            /*             case PTR_LEFT: */
            /*                 new_ref_str += RED(ref->fasta.at(ctg)[beg+ref_ptr++]); */
            /*                 new_alt_str += " "; */
            /*                 new_dels++; */
            /*                 break; */
            /*         } */
            /*     } */

            /*     if (g.print_verbosity >= 2) { */
            /*         printf("\nPTRS: "); */
            /*         for(size_t i = 0; i < cig.size(); i++) { */
            /*             printf("%i ", cig[i]); */
            /*         } */

            /*         printf("\nCIGAR: "); */
            /*         for(size_t i = 0; i < cig.size(); i++) { */
            /*             switch(cig[i]) { */
            /*                 case PTR_DIAG: */
            /*                     if (cig[++i] != PTR_DIAG) */
            /*                         ERROR("cig should be PTR_DIAG"); */
            /*                     printf("M"); break; */
            /*                 case PTR_UP: */
            /*                     printf("I"); break; */
            /*                 case PTR_LEFT: */
            /*                     printf("D"); break; */
            /*             } */
            /*         } */
            /*         printf("\n"); */
            /*     } */

            /*     // print group info, and old/new alignments */
            /*     if (g.print_verbosity >= 1) { */
            /*         if (subs*2 + inss + dels != s) { */
            /*             const char* var_ref; */
            /*             const char* var_alt; */
            /*             printf("\n\n  Group %i: %d variants, %s:%d-%d\n\n" */
            /*                     "    old edit distance: %d (%dS %dI %dD)\n", */
            /*                     int(var_grp), end_idx-beg_idx, ctg.data(), */ 
            /*                     beg, end, subs*2 + inss + dels, subs, inss, dels); */
            /*             for (int var = beg_idx; var < end_idx; var++) { */
            /*                 if (vars.refs[var].size() == 0) var_ref = "-"; */ 
            /*                 else var_ref = vars.refs[var].data(); */
            /*                 if (vars.alts[var].size() == 0) var_alt = "-"; */ 
            /*                 else var_alt = vars.alts[var].data(); */
            /*                 printf("      %s:%i hap%i %s\t%s\t%s\n", */ 
            /*                         ctg.data(), vars.poss[var], h, */ 
            /*                         type_strs[vars.types[var]].data(), */ 
            /*                         var_ref, var_alt); */
            /*             } */
            /*             printf("      REF: %s\n      ALT: %s\n", */ 
            /*                     ref_str.data(), alt_str.data()); */
            /*             printf("\n    new edit distance: %i (%iI %iD)\n", */ 
            /*                     s, new_inss, new_dels); */
            /*             printf("      REF: %s\n      ALT: %s\n", */ 
            /*                     new_ref_str.data(), new_alt_str.data()); */
            /*         } */
            /*     } */

            /*     // update counters */
            /*     old_ed += subs*2 + inss + dels; */
            /*     new_ed += s; */
            /*     if (subs*2 + inss + dels > s) new_ed_groups++; */
            /*     groups++; */

            /* } // group */
        /* } // contig */
    /* } // hap */
    /* INFO("Edit dist reduced in %i of %i groups, from %i to %i.", */ 
    /*         new_ed_groups, groups, old_ed, new_ed); */
    /* return 0; */

/* } */
