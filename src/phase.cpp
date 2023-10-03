#include <vector>
#include <algorithm>

#include "phase.h"
#include "print.h"
#include "globals.h"

void phaseData::write_summary_vcf(std::string out_vcf_fn) {

    // VCF header
    if (g.verbosity >= 1) INFO("  Printing GA4GH-compatible summary VCF to '%s'", out_vcf_fn.data());
    FILE* out_vcf = fopen(out_vcf_fn.data(), "w");
    const std::chrono::time_point now{std::chrono::system_clock::now()};
    time_t tt = std::chrono::system_clock::to_time_t(now);
    tm local_time = *localtime(&tt);
    fprintf(out_vcf, "##fileformat=VCFv4.2\n");
    fprintf(out_vcf, "##fileDate=%04d%02d%02d\n", local_time.tm_year + 1900, 
            local_time.tm_mon + 1, local_time.tm_mday);
    fprintf(out_vcf, "##CL=%s\n", g.cmd.data()+1);
    for (size_t i = 0; i < this->contigs.size(); i++) {
        fprintf(out_vcf, "##contig=<ID=%s,length=%d,ploidy=%d>\n", 
                this->contigs[i].data(), this->lengths[i], this->ploidy[i]);
    }
    fprintf(out_vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GenoType\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Benchmark Decision for call (TP/FP/FN). PP conversion: TP if BC >= 0.5 else FP/FN (hap.py compatibility)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BC,Number=1,Type=String,Description=\"Benchmark Credit for call on the interval [0,1].\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"BenchmarK category (for hap.py compatibility, always genotype match 'gm')\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"variant Quality for ROC creation\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"SuperCluster (index in contig)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=SG,Number=1,Type=Integer,Description=\"Sync Group (for equivalence and credit assignment)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PB,Number=1,Type=Integer,Description=\"Phase Block (index in contig)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Switch (phaseblock state)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PF,Number=1,Type=Integer,Description=\"Phase Flip (supercluster error)\">\n");
    fprintf(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n");

    // write variants
    for (std::string ctg : this->contigs) {
        std::vector< std::vector<int> > ptrs = std::vector< std::vector<int> >(
                CALLSETS, std::vector<int>(HAPS, 0));
        std::vector< std::vector<int> > poss = std::vector< std::vector<int> >(
                CALLSETS, std::vector<int>(HAPS, 0));
        std::vector< std::vector<bool> > next = std::vector< std::vector<bool> >(
                CALLSETS, std::vector<bool>(HAPS, false));
        int sc_idx = 0;
        int ploidy = this->ploidy[std::find(contigs.begin(), contigs.end(), ctg)-contigs.begin()];
        if (this->ctg_phasings[ctg]->sc_phasings.size() == 0) continue;


        // set supercluster flip/swap based on phaseblock and sc phasing
        int phase_block = 0;
        bool phase_switch = this->ctg_phasings[ctg]->pb_phasings[phase_block];
        int phase_sc = this->ctg_phasings[ctg]->sc_phasings[sc_idx];
        bool phase_flip, swap;
        if (phase_switch) {
            if (phase_sc == PHASE_ORIG) { // flip error
                phase_flip = true;
                swap = false;
            } else { // PHASE_SWAP or PHASE_NONE
                phase_flip = false;
                swap = true;
            }
        } else { // no phase switch
            if (phase_sc == PHASE_SWAP) { // flip error
                phase_flip = true;
                swap = true;
            } else { // PHASE_ORIG or PHASE_NONE
                phase_flip = false;
                swap = false;
            }
        }

        auto & vars = this->ctg_phasings[ctg]->ctg_superclusters->ctg_variants;
        while ( ptrs[QUERY][HAP1] < int(vars[QUERY][HAP1]->poss.size()) ||
                ptrs[QUERY][HAP2] < int(vars[QUERY][HAP2]->poss.size()) ||
                ptrs[TRUTH][HAP1] < int(vars[TRUTH][HAP1]->poss.size()) ||
                ptrs[TRUTH][HAP2] < int(vars[TRUTH][HAP2]->poss.size())) {

            // get next positions
            for (int c = 0; c < CALLSETS; c++) {
                for (int h = 0; h < HAPS; h++) {
                    poss[c][h] = ptrs[c][h] < int(vars[c][h]->poss.size()) ? 
                            vars[c][h]->poss[ptrs[c][h]] : 
                            std::numeric_limits<int>::max();
                    if (ptrs[c][h] < int(vars[c][h]->types.size()) &&
                            (vars[c][h]->types[ptrs[c][h]] == TYPE_INS ||
                            vars[c][h]->types[ptrs[c][h]] == TYPE_DEL)) poss[c][h]--;
                }
            }

            // set flags for next haps
            int pos = std::min( std::min(poss[QUERY][HAP1], poss[QUERY][HAP2]),
                    std::min(poss[TRUTH][HAP1], poss[TRUTH][HAP2]));
            for (int c = 0; c < CALLSETS; c++) {
                for (int h = 0; h < HAPS; h++) {
                    next[c][h] = (poss[c][h] == pos);
                }
            }

            // update supercluster and if phase is swapped
            if (pos >= this->ctg_phasings[ctg]->ctg_superclusters->ends[sc_idx]) {
                // update supercluster
                sc_idx++;
                phase_sc = this->ctg_phasings[ctg]->sc_phasings[sc_idx];

                // update phase block
                if (sc_idx >= this->ctg_phasings[ctg]->phase_blocks[phase_block+1])
                    phase_block++;
                phase_switch = this->ctg_phasings[ctg]->pb_phasings[phase_block];

                // update switch/flip status
                if (phase_switch) {
                    if (phase_sc == PHASE_ORIG) { // flip error
                        phase_flip = true;
                        swap = false;
                    } else { // PHASE_SWAP or PHASE_NONE
                        phase_flip = false;
                        swap = true;
                    }
                } else { // no phase switch
                    if (phase_sc == PHASE_SWAP) { // flip error
                        phase_flip = true;
                        swap = true;
                    } else { // PHASE_ORIG or PHASE_NONE
                        phase_flip = false;
                        swap = false;
                    }
                }
            }

            // store if query and truth contain same variant
            std::vector<bool> pair = {false, false};
            for (int h = 0; h < HAPS; h++) {
                pair[h] = next[QUERY][h] && next[TRUTH][swap^h] &&
                        vars[QUERY][h]->refs[ptrs[QUERY][h]] == 
                        vars[TRUTH][swap^h]->refs[ptrs[TRUTH][swap^h]] &&
                        vars[QUERY][h]->alts[ptrs[QUERY][h]] == 
                        vars[TRUTH][swap^h]->alts[ptrs[TRUTH][swap^h]];
            }

            if (next[QUERY][HAP1]) {
                if (next[QUERY][HAP2]) { // two query variants
                    // don't collapse homozygous variants (credit results may differ)
                    for(int h = 0; h < HAPS; h++) {
                        if (pair[h]) { // query matches truth
                            vars[QUERY][h]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY][h]);
                            vars[TRUTH][h^swap]->print_var_sample(out_vcf, ptrs[TRUTH][h^swap], 
                                    (h^swap)?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip);
                            vars[QUERY][h]->print_var_sample(out_vcf, ptrs[QUERY][h], 
                                    h?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip, true);
                            ptrs[QUERY][h]++; ptrs[TRUTH][h^swap]++;
                        } else { // just query
                            vars[QUERY][h]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY][h]);
                            vars[TRUTH][h^swap]->print_var_empty(out_vcf, sc_idx, phase_block);
                            vars[QUERY][h]->print_var_sample(out_vcf, ptrs[QUERY][h], 
                                    h?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip, true);
                            ptrs[QUERY][h]++;
                        }
                    }
                } else { // just query 1
                    if (pair[HAP1]) { // query matches truth
                        vars[QUERY][HAP1]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY][HAP1]);
                        vars[TRUTH][HAP1^swap]->print_var_sample(out_vcf, ptrs[TRUTH][HAP1^swap], 
                                ploidy == 1 ? "1" : (HAP1^swap)?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip);
                        vars[QUERY][HAP1]->print_var_sample(out_vcf, ptrs[QUERY][HAP1], 
                                ploidy == 1 ? "1" : "1|0", sc_idx, phase_block, phase_switch, phase_flip, true);
                        ptrs[QUERY][HAP1]++; ptrs[TRUTH][HAP1^swap]++;
                    } else { // just query
                        vars[QUERY][HAP1]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY][HAP1]);
                        vars[TRUTH][HAP1^swap]->print_var_empty(out_vcf, sc_idx, phase_block);
                        vars[QUERY][HAP1]->print_var_sample(out_vcf, ptrs[QUERY][HAP1], 
                                ploidy == 1 ? "1" : "1|0", sc_idx, phase_block, phase_switch, phase_flip, true);
                        ptrs[QUERY][HAP1]++;
                    }
                }
            } else if (next[QUERY][HAP2]) { // just query 2
                if (pair[HAP2]) { // query matches truth
                    vars[QUERY][HAP2]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY][HAP2]);
                    vars[TRUTH][HAP2^swap]->print_var_sample(out_vcf, ptrs[TRUTH][HAP2^swap], 
                            ploidy == 1 ? "1" : (HAP2^swap)?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip);
                    vars[QUERY][HAP2]->print_var_sample(out_vcf, ptrs[QUERY][HAP2], 
                            ploidy == 1 ? "1" : "0|1", sc_idx, phase_block, phase_switch, phase_flip, true);
                    ptrs[QUERY][HAP2]++; ptrs[TRUTH][HAP2^swap]++;
                } else { // just query
                    vars[QUERY][HAP2]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY][HAP2]);
                    vars[TRUTH][HAP2^swap]->print_var_empty(out_vcf, sc_idx, phase_block);
                    vars[QUERY][HAP2]->print_var_sample(out_vcf, ptrs[QUERY][HAP2], 
                            ploidy == 1 ? "1" : "0|1", sc_idx, phase_block, phase_switch, phase_flip, true);
                    ptrs[QUERY][HAP2]++;
                }
            } else { // no query variants, just truth
                if (next[TRUTH][HAP1] && next[TRUTH][HAP2]) { // two truth variants
                    // don't collapse homozygous variants (credit results may differ)
                    for(int h = 0; h < HAPS; h++) {
                        vars[TRUTH][h]->print_var_info(out_vcf, this->ref, ctg, ptrs[TRUTH][h]);
                        vars[TRUTH][h]->print_var_sample(out_vcf, ptrs[TRUTH][h], 
                                (h^swap)?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip);
                        vars[QUERY][h]->print_var_empty(out_vcf, sc_idx, phase_block, true);
                        ptrs[TRUTH][h]++;
                    }
                } else { // one truth variant
                    if (next[TRUTH][HAP1]) {
                        vars[TRUTH][HAP1]->print_var_info(out_vcf, this->ref, ctg, ptrs[TRUTH][HAP1]);
                        vars[TRUTH][HAP1]->print_var_sample(out_vcf, ptrs[TRUTH][HAP1], 
                                ploidy == 1 ? "1" : swap?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip);
                        vars[QUERY][HAP1]->print_var_empty(out_vcf, sc_idx, phase_block, true);
                        ptrs[TRUTH][HAP1]++;
                    } else if (next[TRUTH][HAP2]) {
                        vars[TRUTH][HAP2]->print_var_info(out_vcf, this->ref, ctg, ptrs[TRUTH][HAP2]);
                        vars[TRUTH][HAP2]->print_var_sample(out_vcf, ptrs[TRUTH][HAP2], 
                                ploidy == 1 ? "1" : swap?"0|1":"1|0", sc_idx, phase_block, phase_switch, phase_flip);
                        vars[QUERY][HAP2]->print_var_empty(out_vcf, sc_idx, phase_block, true);
                        ptrs[TRUTH][HAP2]++;
                    } else {
                        ERROR("No variants are selected next.");
                    }
                }
            }
        }
    }
    fclose(out_vcf);
}


/*******************************************************************************/


phaseData::phaseData(std::shared_ptr<superclusterData> clusterdata_ptr)
{
    // copy contigs and reference
    for (int i = 0; i < int(clusterdata_ptr->contigs.size()); i++) {
        std::string ctg = clusterdata_ptr->contigs[i];
        this->contigs.push_back(ctg);
        this->lengths.push_back(clusterdata_ptr->lengths[i]);
        this->ploidy.push_back(clusterdata_ptr->ploidy[i]);
        this->ctg_phasings[ctg] = std::shared_ptr<ctgPhasings>(new ctgPhasings());
    }
    this->ref = clusterdata_ptr->ref;

    // for each contig
    for (auto ctg : clusterdata_ptr->contigs) {

        // add pointer to clusters
        this->ctg_phasings[ctg]->ctg_superclusters = clusterdata_ptr->ctg_superclusters[ctg];

        // add ctg_phasings
        for(int i = 0; i < clusterdata_ptr->ctg_superclusters[ctg]->n; i++) {
            this->ctg_phasings[ctg]->sc_phasings.push_back(
                    clusterdata_ptr->ctg_superclusters[ctg]->phase[i]);
            this->ctg_phasings[ctg]->n++;
        }
    }
    this->phase();
}

void phaseData::phase()
{
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[7/8] Phasing superclusters%s",
            COLOR_PURPLE, COLOR_WHITE);
    // phase each contig separately
    for (auto ctg : this->contigs) {
        std::vector< std::vector<int> > 
            mat(2, std::vector<int>(this->ctg_phasings[ctg]->n+1));
        std::vector< std::vector<int> > 
            ptrs(2, std::vector<int>(this->ctg_phasings[ctg]->n));

        // forward pass
        for (int i = 0; i < this->ctg_phasings[ctg]->n; i++) {

            // determine costs (penalized if this phasing deemed incorrect)
            std::vector<int> costs = {0, 0};
            switch (this->ctg_phasings[ctg]->sc_phasings[i]) {
                case PHASE_NONE: 
                    costs[PHASE_PTR_KEEP] = 0; 
                    costs[PHASE_PTR_SWAP] = 0; 
                    break;
                case PHASE_ORIG: 
                    costs[PHASE_PTR_KEEP] = 0; 
                    costs[PHASE_PTR_SWAP] = 1; 
                     break;
                case PHASE_SWAP:
                    costs[PHASE_PTR_KEEP] = 1; 
                    costs[PHASE_PTR_SWAP] = 0; 
                     break;
                default:
                    ERROR("Unexpected phase (%d)", this->ctg_phasings[ctg]->sc_phasings[i]);
                    break;
            }

            for (int phase = 0; phase < 2; phase++) {
                if (mat[phase][i] + costs[phase] <= mat[phase^1][i] + 1) {
                    mat[phase][i+1] = mat[phase][i] + costs[phase];
                    ptrs[phase][i] = PHASE_PTR_KEEP;
                }
                else {
                    mat[phase][i+1] = mat[phase^1][i] + 1;
                    ptrs[phase][i] = PHASE_PTR_SWAP;
                }
            }
        }

        // backwards pass
        if (this->ctg_phasings[ctg]->n > 0) { // skip empty contigs
            // determine starting phase
            int phase = 0;
            if (mat[PHASE_SWAP][this->ctg_phasings[ctg]->n] <
                    mat[PHASE_ORIG][this->ctg_phasings[ctg]->n])
                phase = 1;

            int i = this->ctg_phasings[ctg]->n-1;
            this->ctg_phasings[ctg]->phase_blocks.push_back(i+1);
            while (i > 0) {
                if (ptrs[phase][i] == PHASE_PTR_SWAP) { // new phase block
                    this->ctg_phasings[ctg]->nswitches++;
                    this->ctg_phasings[ctg]->phase_blocks.push_back(i+1);
                    this->ctg_phasings[ctg]->pb_phasings.push_back(phase);
                    phase ^= 1;
                } else if (ptrs[phase][i] == PHASE_PTR_KEEP) { // within phase block
                    if (this->ctg_phasings[ctg]->sc_phasings[i] != PHASE_NONE &&
                            this->ctg_phasings[ctg]->sc_phasings[i] != phase)
                        this->ctg_phasings[ctg]->nerrors++;
                }
                i--;
            }
            this->ctg_phasings[ctg]->phase_blocks.push_back(0);
            std::reverse(this->ctg_phasings[ctg]->phase_blocks.begin(),
                    this->ctg_phasings[ctg]->phase_blocks.end());
            this->ctg_phasings[ctg]->pb_phasings.push_back(phase);
            std::reverse(this->ctg_phasings[ctg]->pb_phasings.begin(),
                    this->ctg_phasings[ctg]->pb_phasings.end());
        }
    }
    

    // print
    int switch_errors = 0;
    int phase_errors = 0;
    for (auto ctg : this->contigs) {

        if (g.verbosity >= 2) {
            int s = 0;
            std::vector<std::string> colors {"\033[34m", "\033[32m"};
            int color_idx = 0;
            for (int i = 0; i < this->ctg_phasings[ctg]->n; i++) {
                if (i == this->ctg_phasings[ctg]->phase_blocks[s]) {
                    printf("%s", colors[color_idx].data());
                    color_idx ^= 1;
                    s++;
                }
                int ctg_phasing = this->ctg_phasings[ctg]->sc_phasings[i];
                printf("%s", phase_strs[ctg_phasing].data());
            }
            printf("\033[0m\n");
        }
        switch_errors += this->ctg_phasings[ctg]->nswitches;
        phase_errors += this->ctg_phasings[ctg]->nerrors;
    }
    if (g.verbosity >= 1) INFO("  Total switch errors: %d", switch_errors);
    if (g.verbosity >= 1) INFO("  Total phase errors: %d", phase_errors);
}
