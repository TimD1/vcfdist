#include <vector>
#include <algorithm>
#include <stdexcept>

#include "phase.h"
#include "print.h"
#include "globals.h"


void phaseblockData::write_summary_vcf(std::string out_vcf_fn) {

    // VCF header
    if (g.verbosity >= 1) INFO("  Writing summary VCF to '%s'", out_vcf_fn.data());
    FILE* out_vcf = fopen(out_vcf_fn.data(), "w");
    const std::chrono::time_point<std::chrono::system_clock> now{std::chrono::system_clock::now()};
    time_t tt = std::chrono::system_clock::to_time_t(now);
    tm local_time = *localtime(&tt);
    fprintf(out_vcf, "##fileformat=VCFv4.2\n");
    fprintf(out_vcf, "##fileDate=%04d%02d%02d\n", local_time.tm_year + 1900, 
            local_time.tm_mon + 1, local_time.tm_mday);
    fprintf(out_vcf, "##CL=%s\n", g.cmd.data());
    for (size_t i = 0; i < this->contigs.size(); i++) {
        fprintf(out_vcf, "##contig=<ID=%s,length=%d,ploidy=%d>\n", 
                this->contigs[i].data(), this->lengths[i], this->ploidy[i]);
    }
    fprintf(out_vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GenoType\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Benchmark Decision for call (TP/FP/FN).\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BC,Number=1,Type=Float,Description=\"Benchmark Credit (on the interval [0,1], based on sync group edit distance)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference edit Distance from truth within current sync group\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=QD,Number=1,Type=Integer,Description=\"Query edit Distance from truth within current sync group\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"BenchmarK category ('gm' if credit == 1, 'lm' if credit > 0, else '.')\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"variant Quality\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"SuperCluster (index in contig)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=SG,Number=1,Type=Integer,Description=\"Sync Group (index in supercluster, for credit assignment)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set identifier (input, per-variant)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PB,Number=1,Type=Integer,Description=\"Phase Block (output, per-supercluster, index in contig)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BS,Number=1,Type=Integer,Description=\"Block Phase: 0 = PHASE_KEEP, 1 = PHASE_SWAP)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=VP,Number=1,Type=Integer,Description=\"Variant Phase: 0 = PHASE_ORIG, 1 = PHASE_SWAP, . = PHASE_NONE)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=FE,Number=1,Type=Integer,Description=\"Flip Error (a per-supercluster error)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GE,Number=1,Type=String,Description=\"Genotype Error ('+' if 0/1 truth -> 1/1 query, '-' if 1/1 truth -> 0/1 query, '.' otherwise)\">\n");
    fprintf(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n");

    // write variants
    for (std::string ctg : this->contigs) {
        std::vector<int> ptrs = std::vector<int>(CALLSETS, 0);
        std::vector<int> poss = std::vector<int>(CALLSETS, 0);
        std::vector<int> next = std::vector<int>(CALLSETS, 0);
        int sc_idx = 0;
        int ploidy = this->ploidy[std::find(contigs.begin(), contigs.end(), ctg)-contigs.begin()];
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
        auto & vars = ctg_pbs->ctg_superclusters->callset_vars;
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];
        std::shared_ptr<ctgVariants> tvars = ctg_pbs->ctg_superclusters->callset_vars[TRUTH];
        if (qvars->n == 0) continue;

        // set supercluster flip/swap based on phaseblock and sc phasing
        int phase_block = 0;
        bool block_state = qvars->pb_phases[ptrs[QUERY]];
        int phase = qvars->phases[ptrs[QUERY]];
        bool flip_error;
        if (block_state == PHASE_SWAP) {
            if (phase == PHASE_ORIG) { // flip error
                flip_error = true;
            } else { // PHASE_SWAP or PHASE_NONE
                flip_error = false;
            }
        } else { // block_state == PHASE_ORIG
            if (phase == PHASE_SWAP) { // flip error
                flip_error = true;
            } else { // PHASE_ORIG or PHASE_NONE
                flip_error = false;
            }
        }

        while ( ptrs[QUERY] < qvars->n || ptrs[TRUTH] < tvars->n) {

            // get next positions
            for (int c = 0; c < CALLSETS; c++) {
                poss[c] = ptrs[c] < int(vars[c]->poss.size()) ? 
                        vars[c]->poss[ptrs[c]] : 
                        std::numeric_limits<int>::max();
                if (ptrs[c] < int(vars[c]->types.size()) &&
                        (vars[c]->types[ptrs[c]] == TYPE_INS ||
                        vars[c]->types[ptrs[c]] == TYPE_DEL)) poss[c]--;
            }

            // set flags for next haps
            int pos = std::min(poss[QUERY], poss[TRUTH]);
            for (int c = 0; c < CALLSETS; c++) {
                next[c] = (poss[c] == pos);
            }

            // update phasing
            if (ptrs[QUERY] < qvars->n) {
                phase = qvars->phases[ptrs[QUERY]];
                block_state = qvars->pb_phases[ptrs[QUERY]];

                // update switch/flip status
                if (block_state == PHASE_SWAP) {
                    if (phase == PHASE_ORIG) {
                        flip_error = true;
                    } else { // PHASE_SWAP or PHASE_NONE
                        flip_error = false;
                    }
                } else { // block_state == PHASE_ORIG
                    if (phase == PHASE_SWAP) {
                        flip_error = true;
                    } else { // PHASE_ORIG or PHASE_NONE
                        flip_error = false;
                    }
                }
            }

            // update supercluster and phase block
            if (pos >= ctg_scs->ends[sc_idx]) {
                sc_idx++;
                if (ptrs[QUERY] >= ctg_pbs->phase_blocks[phase_block+1])
                    phase_block++;
            }

            /* if (next[QUERY] && ptrs[QUERY] < qvars->n) { */
            /*     fprintf(out_vcf, "orig_gt: %s\tcalc_gt: %s\tcredit: %.2f|%.2f\tref_dist: %d|%d\n", */ 
            /*             gt_strs[vars[QUERY]->orig_gts[ptrs[QUERY]]].data(), */
            /*             gt_strs[vars[QUERY]->calc_gts[ptrs[QUERY]]].data(), */
            /*             vars[QUERY]->credit[HAP1][ptrs[QUERY]], */
            /*             vars[QUERY]->credit[HAP2][ptrs[QUERY]], */
            /*             vars[QUERY]->ref_ed[HAP1][ptrs[QUERY]], */
            /*             vars[QUERY]->ref_ed[HAP2][ptrs[QUERY]] */
            /*             ); */
            /* } */
            /* if (next[TRUTH] && ptrs[TRUTH] < tvars->n) { */
            /*     fprintf(out_vcf, "orig_gt: %s\n", */ 
            /*             gt_strs[vars[TRUTH]->orig_gts[ptrs[TRUTH]]].data()); */
            /* } */

            if (next[QUERY]) {
                if (next[TRUTH]) {
                    if (vars[QUERY]->refs[ptrs[QUERY]] == vars[TRUTH]->refs[ptrs[TRUTH]] &&
                        vars[QUERY]->alts[ptrs[QUERY]] == vars[TRUTH]->alts[ptrs[TRUTH]]) { // query matches truth
                        // print data for each haplotype
                        for (int qhi = 0; qhi < HAPS; qhi++) {
                            bool swap = vars[QUERY]->calcgt_is_swapped(ptrs[QUERY]);
                            int thi = qhi ^ swap ^ block_state ^ flip_error;
                            if (vars[QUERY]->var_on_hap(ptrs[QUERY], qhi, true) || 
                                    vars[TRUTH]->var_on_hap(ptrs[TRUTH], thi)) {
                                vars[QUERY]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY]);
                                if (vars[TRUTH]->var_on_hap(ptrs[TRUTH], thi)) { // print truth
                                    vars[TRUTH]->print_var_sample(out_vcf, ptrs[TRUTH], thi,
                                        ploidy == 1 ? "1" : (thi ? "0|1" : "1|0"), 
                                        sc_idx, phase_block, block_state, flip_error);
                                } else {
                                    vars[TRUTH]->print_var_empty(out_vcf, sc_idx, phase_block);
                                }
                                if (vars[QUERY]->var_on_hap(ptrs[QUERY], qhi, true)) { // print query
                                    vars[QUERY]->print_var_sample(out_vcf, ptrs[QUERY], qhi,
                                        ploidy == 1 ? "1" : ((qhi ^ swap) ? "0|1" : "1|0"),
                                        sc_idx, phase_block, block_state, flip_error, true);
                                } else {
                                    vars[QUERY]->print_var_empty(out_vcf, sc_idx, phase_block, true);
                                }
                            }
                        }
                        ptrs[QUERY]++; ptrs[TRUTH]++;
                    } else { // positional tie, diff vars, just print query
                        for (int qhi = 0; qhi < HAPS; qhi++) {
                            bool swap = vars[QUERY]->calcgt_is_swapped(ptrs[QUERY]);
                            if (vars[QUERY]->var_on_hap(ptrs[QUERY], qhi, true)) {
                                vars[QUERY]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY]);
                                vars[TRUTH]->print_var_empty(out_vcf, sc_idx, phase_block);
                                vars[QUERY]->print_var_sample(out_vcf, ptrs[QUERY], qhi,
                                        ploidy == 1 ? "1" : ((qhi ^ swap) ? "0|1" : "1|0"),
                                        sc_idx, phase_block, block_state, flip_error, true);
                            }
                        }
                        ptrs[QUERY]++;
                    }
                } else { // query is next
                    for (int qhi = 0; qhi < HAPS; qhi++) {
                        bool swap = vars[QUERY]->calcgt_is_swapped(ptrs[QUERY]);
                        if (vars[QUERY]->var_on_hap(ptrs[QUERY], qhi, true)) {
                            vars[QUERY]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY]);
                            vars[TRUTH]->print_var_empty(out_vcf, sc_idx, phase_block);
                            vars[QUERY]->print_var_sample(out_vcf, ptrs[QUERY], qhi,
                                    ploidy == 1 ? "1" : ((qhi ^ swap) ? "0|1" : "1|0"),
                                    sc_idx, phase_block, block_state, flip_error, true);
                        }
                    }
                    ptrs[QUERY]++;
                }
            } else if (next[TRUTH]) {
                for (int qhi = 0; qhi < HAPS; qhi++) {
                    int thi = qhi ^ block_state ^ flip_error;
                    if (vars[TRUTH]->var_on_hap(ptrs[TRUTH], thi)) {
                        vars[TRUTH]->print_var_info(out_vcf, this->ref, ctg, ptrs[TRUTH]);
                        vars[TRUTH]->print_var_sample(out_vcf, ptrs[TRUTH], thi,
                                ploidy == 1 ? "1" : (thi ? "0|1" : "1|0"), 
                                sc_idx, phase_block, block_state, flip_error);
                        vars[QUERY]->print_var_empty(out_vcf, sc_idx, phase_block, true);
                    }
                }
                ptrs[TRUTH]++;
            } else {
                ERROR("No variants are selected next.");
            }
        }
    }
    fclose(out_vcf);
}


/*******************************************************************************/


phaseblockData::phaseblockData(std::shared_ptr<superclusterData> clusterdata_ptr)
{
    // copy contigs and reference
    for (int i = 0; i < int(clusterdata_ptr->contigs.size()); i++) {
        std::string ctg = clusterdata_ptr->contigs[i];
        this->contigs.push_back(ctg);
        this->lengths.push_back(clusterdata_ptr->lengths[i]);
        this->ploidy.push_back(clusterdata_ptr->ploidy[i]);
        this->phase_blocks[ctg] = std::shared_ptr<ctgPhaseblocks>(new ctgPhaseblocks());
    }
    this->ref = clusterdata_ptr->ref;

    // add pointers to clusters
    for (const std::string & ctg : this->contigs) {
        this->phase_blocks[ctg]->ctg_superclusters = clusterdata_ptr->superclusters[ctg];
    }

    // add phase blocks based on phase sets for each contig
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];
        int curr_phase_set = -1;
        for (int qvar_idx = 0; qvar_idx < qvars->n; qvar_idx++) {
            if (qvars->phase_sets[qvar_idx] != curr_phase_set) {
                ctg_pbs->phase_blocks.push_back(qvar_idx);
                ctg_pbs->n++;
                curr_phase_set = qvars->phase_sets[qvar_idx];
            }
        }
        ctg_pbs->phase_blocks.push_back(qvars->n);
    }
    
    // calculate phasings, flip, and switch errors
    this->phase();

    // fix genotypes, with ties defaulting to current phase
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgSuperclusters> scs = clusterdata_ptr->superclusters[ctg];
        for (int sc_idx = 0; sc_idx < scs->n; sc_idx++) {
            fix_genotype_allele_counts(scs, sc_idx);
        }
    }
}


/*******************************************************************************/


void phaseblockData::phase()
{
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[6/7] Phasing superclusters%s",
            COLOR_PURPLE, COLOR_WHITE);

    // phase each contig separately
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];
        std::vector< std::vector<int> > mat(2, std::vector<int>(qvars->n+1));
        std::vector< std::vector<int> > ptrs(2, std::vector<int>(qvars->n+1));

        // calculate phasings for each variant
        for (int i = 0; i < qvars->n; i++) {
            if ((qvars->orig_gts[i] == GT_ALT1_REF && qvars->calc_gts[i] == GT_ALT1_REF) || // same
                    (qvars->orig_gts[i] == GT_REF_ALT1 && qvars->calc_gts[i] == GT_REF_ALT1)) {
                qvars->phases[i] = PHASE_ORIG;
            } else if ((qvars->orig_gts[i] == GT_ALT1_REF && qvars->calc_gts[i] == GT_REF_ALT1) || // diff
                    (qvars->orig_gts[i] == GT_REF_ALT1 && qvars->calc_gts[i] == GT_ALT1_REF)) {
                qvars->phases[i] = PHASE_SWAP;
            } else {
                qvars->phases[i] = PHASE_NONE;
            }
        }

        // forward pass
        for (int i = 0; i < qvars->n; i++) {

            // determine costs (penalized if this phasing deemed incorrect)
            std::vector<int> costs = {0, 0};
            switch (qvars->phases[i]) {
                case PHASE_ORIG:
                    costs[PHASE_ORIG] = 0; 
                    costs[PHASE_SWAP] = 1; 
                    break;
                case PHASE_SWAP:
                    costs[PHASE_ORIG] = 1; 
                    costs[PHASE_SWAP] = 0; 
                    break;
                case PHASE_NONE:
                    costs[PHASE_ORIG] = 0; 
                    costs[PHASE_SWAP] = 0; 
                    break;
                default:
                    ERROR("Unexpected phase (%d)", qvars->phases[i]);
                    break;
            }

            // mat[i] ptrs[i] refer to points directly before cluster i, and if switch occurs it
            // is between cluster i and i+1

            // no cost for phase switches if on border of phase set
            int cost_swap = 1;
            if (i < qvars->n-1 && qvars->phase_sets[i] != qvars->phase_sets[i+1]) {
                cost_swap = 0;
            }

            for (int phase = 0; phase < 2; phase++) {
                if (mat[phase][i] + costs[phase] < mat[phase^1][i] + costs[phase^1] + cost_swap) {
                    mat[phase][i+1] = mat[phase][i] + costs[phase];
                    ptrs[phase][i+1] = PHASE_PTR_KEEP;
                }
                else {
                    mat[phase][i+1] = mat[phase^1][i] + costs[phase^1] + cost_swap;
                    ptrs[phase][i+1] = PHASE_PTR_SWAP;
                }
            }
        }

        // backwards pass
        if (qvars->n > 0) { // skip empty contigs

            // determine starting phase
            int phase = PHASE_ORIG; // 0
            if (mat[PHASE_SWAP][qvars->n] < mat[PHASE_ORIG][qvars->n])
                phase = PHASE_SWAP; // 1

            int i = qvars->n;
            while (i > 0) {
                if (ptrs[phase][i] == PHASE_PTR_SWAP) {

                    // not a switch error if between phase sets
                    if (qvars->phase_sets[i] == qvars->phase_sets[i-1]) {
                        ctg_pbs->switches.push_back(i);
                        ctg_pbs->nswitches++;
                    }

                    phase ^= 1;
                } else if (ptrs[phase][i] == PHASE_PTR_KEEP) { // within phase block
                    if (qvars->phases[i-1] != PHASE_NONE && qvars->phases[i-1] != phase) {
                        ctg_pbs->flips.push_back(i-1);
                        ctg_pbs->nflips++;
                    }
                }
                i--;
                qvars->pb_phases[i] = phase;
            }
            std::reverse(ctg_pbs->flips.begin(), ctg_pbs->flips.end());
            std::reverse(ctg_pbs->switches.begin(), ctg_pbs->switches.end());
        }
    }
    

    // print
    int switch_errors = 0;
    int flip_errors = 0;
    int phase_blocks = 0;
    int variants = 0;
    if (g.verbosity >= 1) INFO("  Contigs:");
    int id = 0;
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];

        // print errors per contig
        if (g.verbosity >= 1) {
            INFO("    [%2d] %s: %d switch errors, %d flip errors, %d phase blocks", id, ctg.data(), 
                    ctg_pbs->nswitches, ctg_pbs->nflips, ctg_pbs->n);
        }
        variants += qvars->n;
        switch_errors += ctg_pbs->nswitches;
        flip_errors += ctg_pbs->nflips;
        phase_blocks += ctg_pbs->n;
        id++;
    }
    int ng50 = this->calculate_ng50(false, false);
    int s_ngc50 = this->calculate_ng50(true, false);
    int sf_ngc50 = this->calculate_ng50(true, true);

    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("             Total  phase blocks: %d", phase_blocks);
    if (g.verbosity >= 1) INFO("             Total switch errors: %d", switch_errors);
    if (g.verbosity >= 1) INFO("             Total   flip errors: %d", flip_errors);
    if (g.verbosity >= 1 && variants) 
        INFO("               Switch error rate: %.6f%%", 100*switch_errors/float(variants));
    if (g.verbosity >= 1 && variants) 
        INFO("                 Flip error rate: %.6f%%", 100*flip_errors/float(variants));
    if (g.verbosity >= 1) INFO("               Phase block  NG50: %d", ng50);
    if (g.verbosity >= 1) INFO("  (switch)     Phase block NGC50: %d", s_ngc50);
    if (g.verbosity >= 1) INFO("  (switchflip) Phase block NGC50: %d", sf_ngc50);
    this->write_phasing_summary(phase_blocks, switch_errors, flip_errors, variants, ng50, s_ngc50, sf_ngc50);
}


/*******************************************************************************/


void phaseblockData::write_switchflips() {

    std::string out_sf_fn = g.out_prefix + "switchflips.tsv";
    FILE* out_sf = 0;
    if (g.verbosity >= 1) INFO("  Writing switchflip results to '%s'", out_sf_fn.data());
    out_sf = fopen(out_sf_fn.data(), "w");
    fprintf(out_sf, "CONTIG\tSTART\tSTOP\tSWITCH_TYPE\tVARIANT\tPHASE_BLOCK\n");

    // get sizes of each correct phase block (split on flips, not just switch)
    for (std::string ctg: this->contigs) {
        
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];

        int switch_idx = 0;
        int flip_idx = 0;
        int pb_idx = 1;
        if (qvars->n == 0) continue;

        // init start of this correct block
        int vi = 0;
        int next_vi = qvars->n; // default to last
        int beg = 0; int end = 0;
        int type = SWITCHTYPE_NONE;

        while (true) {

            // get type and supercluster of next switch/flip/phaseset
            if (pb_idx < ctg_pbs->n && ctg_pbs->phase_blocks[pb_idx] <= next_vi) {
                type = SWITCHTYPE_SWITCH;
                next_vi = ctg_pbs->phase_blocks[pb_idx];
            }
            if (switch_idx < ctg_pbs->nswitches && ctg_pbs->switches[switch_idx] <= next_vi) {
                type = SWITCHTYPE_SWITCH_ERR;
                next_vi = ctg_pbs->switches[switch_idx];
            }
            // check flip last because it will cause a second switch
            if (flip_idx < ctg_pbs->nflips && ctg_pbs->flips[flip_idx] <= next_vi) {
                if (type == SWITCHTYPE_SWITCH)
                    type = SWITCHTYPE_SWITCH_AND_FLIP;
                else
                    type = SWITCHTYPE_FLIP;
                next_vi = ctg_pbs->flips[flip_idx];
            }
            if (type == SWITCHTYPE_NONE) { // all out-of-bounds
                break;
            }
            if (next_vi < vi) ERROR("Next variant (%d) is not after current variant (%d) in write_switchflips()", next_vi, vi);


            // get block(s)
            if (type == SWITCHTYPE_FLIP || type == SWITCHTYPE_SWITCH_AND_FLIP) {
                // switch could have occurred anywhere after last phased supercluster
                int left = next_vi-1;
                while (left > 0 && qvars->phases[left] == PHASE_NONE)
                    left--;
                if (left >= 0) {
                    beg = qvars->poss[left] + qvars->rlens[left];
                    end = qvars->poss[next_vi];
                    fprintf(out_sf, "%s\t%d\t%d\t%s\t%d\t%d\n", ctg.data(), beg, end, 
                            switch_strs[SWITCHTYPE_FLIP_BEG].data(), next_vi, pb_idx-1);
                }

                // switch could have occurred anywhere before next phased supercluster
                int right = next_vi+1;
                while (right < qvars->n-1 && qvars->phases[right] == PHASE_NONE)
                    right++;
                if (right < qvars->n) {
                    beg = qvars->poss[next_vi] + qvars->rlens[next_vi];
                    end = qvars->poss[right];
                    fprintf(out_sf, "%s\t%d\t%d\t%s\t%d\t%d\n", ctg.data(), beg, end, 
                            switch_strs[SWITCHTYPE_FLIP_END].data(), next_vi, pb_idx-1);
                }
                flip_idx++;
                if (type == SWITCHTYPE_SWITCH_AND_FLIP) pb_idx++;

            } else if (type == SWITCHTYPE_SWITCH) {
                // end of phase block, don't print anything since not an error
                pb_idx++;

            } else if (type == SWITCHTYPE_SWITCH_ERR) {
                // expand left/right from in between these clusters
                int left = next_vi-1;
                while (left > 0 && qvars->phases[left] == PHASE_NONE)
                    left--;
                int right = next_vi;
                while (right < qvars->n-1 && qvars->phases[right] == PHASE_NONE)
                    right++;
                if (left >= 0 && right < qvars->n) {
                    beg = qvars->poss[left] + qvars->rlens[left];
                    end = qvars->poss[right];
                    fprintf(out_sf, "%s\t%d\t%d\t%s\t%d\t%d\n", ctg.data(), beg, end, 
                            switch_strs[SWITCHTYPE_SWITCH_ERR].data(), next_vi, pb_idx-1);
                }
                switch_idx++;
            }

            vi = next_vi;
            next_vi = qvars->n; // reset to end
            type = SWITCHTYPE_NONE;
        }
    }
    fclose(out_sf);
}


/*******************************************************************************/


void phaseblockData::write_phasing_summary(int phase_blocks, int switch_errors,
        int flip_errors, int variants, int ng50, int s_ngc50, int sf_ngc50) {
    std::string out_phasing_summary_fn = g.out_prefix + "phasing-summary.tsv";
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Writing phasing summary to '%s'", 
            out_phasing_summary_fn.data());
    FILE* out_phasing_summary = fopen(out_phasing_summary_fn.data(), "w");
    fprintf(out_phasing_summary,
            "PHASE_BLOCKS\tSWITCH_ERRORS\tFLIP_ERRORS\tSWITCH_ERROR_RATE\tFLIP_ERROR_RATE\t"
            "NG_50\tSWITCH_NGC50\tSWITCHFLIP_NGC50\n");
    fprintf(out_phasing_summary, "%d\t%d\t%d\t%.6f%%\t%.6f%%\t%d\t%d\t%d", phase_blocks, 
            switch_errors, flip_errors, 100*switch_errors/float(variants),
            100*flip_errors/float(variants), ng50, s_ngc50, sf_ngc50);
    fclose(out_phasing_summary);
}


/*******************************************************************************/


int phaseblockData::calculate_ng50(bool break_on_switch, bool break_on_flip) {

    // get total bases in genome
    size_t total_bases = 0;
    for (size_t i = 0; i < this->contigs.size(); i++) {
        total_bases += lengths[i];
    }

    // get sizes of each correct phase block (split on flips, not just switch)
    std::vector<int> correct_blocks;
    for (std::string ctg: this->contigs) {
        
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];

        int switch_idx = 0;
        int flip_idx = 0;
        int pb_idx = 1;
        if (qvars->n == 0) continue;

        // init start of this correct block
        int vi = 0;
        int next_vi = qvars->n; // default to last
        int beg = qvars->poss[0];
        int end = 0;
        int type = SWITCHTYPE_NONE;

        while (true) {

            // get type and supercluster of next switch/flip/phaseset
            if (pb_idx < ctg_pbs->n && ctg_pbs->phase_blocks[pb_idx] <= next_vi) {
                type = SWITCHTYPE_SWITCH;
                next_vi = ctg_pbs->phase_blocks[pb_idx];
            }
            if (break_on_switch && switch_idx < ctg_pbs->nswitches && ctg_pbs->switches[switch_idx] <= next_vi) {
                type = SWITCHTYPE_SWITCH_ERR;
                next_vi = ctg_pbs->switches[switch_idx];
            }
            // check flip last (takes preference due to <=) because it can cause two breaks (before/after)
            // NOTE: is is possible for one supercluster to have both a switch (new PS) and flip
            if (break_on_flip && flip_idx < ctg_pbs->nflips && ctg_pbs->flips[flip_idx] <= next_vi) {
                if (type == SWITCHTYPE_SWITCH)
                    type = SWITCHTYPE_SWITCH_AND_FLIP;
                else
                    type = SWITCHTYPE_FLIP;
                next_vi = ctg_pbs->flips[flip_idx];
            }
            if (type == SWITCHTYPE_NONE) { // all out-of-bounds
                break;
            }
            if (next_vi < vi) ERROR("Next variant (%d) is not after current variant (%d) in calc_ng50()", next_vi, vi);


            // get block(s)
            if (type == SWITCHTYPE_FLIP || type == SWITCHTYPE_SWITCH_AND_FLIP) {
                end = qvars->poss[next_vi-1] + qvars->rlens[next_vi-1];
                correct_blocks.push_back(end-beg);
                beg = qvars->poss[next_vi];

                end = qvars->poss[next_vi] + qvars->rlens[next_vi];
                correct_blocks.push_back(end-beg);
                beg = qvars->poss[next_vi+1];
                flip_idx++;
                if (type == SWITCHTYPE_SWITCH_AND_FLIP) pb_idx++;
            } else if (type == SWITCHTYPE_SWITCH) {
                end = qvars->poss[next_vi-1] + qvars->rlens[next_vi-1];
                correct_blocks.push_back(end-beg);
                beg = qvars->poss[next_vi];
                pb_idx++;
            } else if (type == SWITCHTYPE_SWITCH_ERR) {
                end = qvars->poss[next_vi-1] + qvars->rlens[next_vi-1];
                correct_blocks.push_back(end-beg);
                beg = qvars->poss[next_vi];
                switch_idx++;
            }

            vi = next_vi;
            next_vi = qvars->n; // reset to end
            type = SWITCHTYPE_NONE;
        }
        end = qvars->poss[qvars->n-1] + qvars->rlens[qvars->n-1];
        correct_blocks.push_back(end-beg);
    }

    // return NGC50
    size_t total_correct = 0;
    std::sort(correct_blocks.begin(), correct_blocks.end(), std::greater<>());
    for (size_t i = 0; i < correct_blocks.size(); i++) {
        total_correct += correct_blocks[i];
        if (total_correct >= total_bases / 2)
            return correct_blocks[i];
    }
    return 0;
}
