#include <vector>
#include <algorithm>

#include "phase.h"
#include "print.h"
#include "globals.h"

phaseData::phaseData(std::shared_ptr<superclusterData> clusterdata_ptr)
{
    // copy contigs
    for (auto ctg : clusterdata_ptr->contigs) {
        this->contigs.push_back(ctg);
        ctg_phasings[ctg] = ctgPhasings();
    }

    // for each contig
    for (auto ctg : clusterdata_ptr->contigs) {

        // add pointer to clusters
        this->ctg_phasings[ctg].ctg_superclusters = clusterdata_ptr->ctg_superclusters[ctg];

        // add ctg_phasings
        for(int i = 0; i < clusterdata_ptr->ctg_superclusters[ctg]->n; i++) {
            this->ctg_phasings[ctg].phasings.push_back(
                    clusterdata_ptr->ctg_superclusters[ctg]->phase[i]);
            this->ctg_phasings[ctg].n++;
        }
    }
    this->phase();
}

void phaseData::phase()
{
    INFO(" ");
    INFO("Phasing superclusters");
    // phase each contig separately
    for (auto ctg : this->contigs) {
        std::vector< std::vector<int> > 
            mat(2, std::vector<int>(this->ctg_phasings[ctg].n+1));
        std::vector< std::vector<int> > 
            ptrs(2, std::vector<int>(this->ctg_phasings[ctg].n));
        if (this->ctg_phasings[ctg].n)
            INFO("  Contig '%s'", ctg.data());

        // forward pass
        for (int i = 0; i < this->ctg_phasings[ctg].n; i++) {

            // determine costs (penalized if this phasing deemed incorrect)
            std::vector<int> costs = {0, 0};
            switch (this->ctg_phasings[ctg].phasings[i]) {
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
                    ERROR("Unexpected phase (%d)", this->ctg_phasings[ctg].phasings[i]);
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
        if (this->ctg_phasings[ctg].n > 0) { // skip empty contigs
            int i = this->ctg_phasings[ctg].n-1;
            int phase = 0;
            this->ctg_phasings[ctg].phase_blocks.push_back(i+1);
            while (i > 0) {
                if (ptrs[phase][i] == PHASE_PTR_SWAP) {
                    phase ^= 1;
                    this->ctg_phasings[ctg].nswitches++;
                    this->ctg_phasings[ctg].phase_blocks.push_back(i+1);
                }
                i--;
            }
            this->ctg_phasings[ctg].phase_blocks.push_back(0);
            std::reverse(this->ctg_phasings[ctg].phase_blocks.begin(),
                    this->ctg_phasings[ctg].phase_blocks.end());
        }
        if (this->ctg_phasings[ctg].n)
            INFO("    %d switch errors", this->ctg_phasings[ctg].nswitches);
    }
    

    // print
    int switch_errors = 0;
    for (auto ctg : this->contigs) {

        if (g.print_verbosity >= 2) {
            int s = 0;
            std::vector<std::string> colors {"\033[34m", "\033[32m"};
            int color_idx = 0;
            for (int i = 0; i < this->ctg_phasings[ctg].n; i++) {
                if (i == this->ctg_phasings[ctg].phase_blocks[s]) {
                    printf("%s", colors[color_idx].data());
                    color_idx ^= 1;
                    s++;
                }
                int ctg_phasing = this->ctg_phasings[ctg].phasings[i];
                printf("%s", phase_strs[ctg_phasing].data());
            }
            printf("\033[0m\n");
        }

        if (g.print_verbosity >= 1) {
            printf("Contig '%s' phase block cluster indices: ", ctg.data());
            for(int i = 0; i < this->ctg_phasings[ctg].nswitches; i++)
                printf("%d ", this->ctg_phasings[ctg].phase_blocks[i]);
            printf("\n");
        }
        switch_errors += this->ctg_phasings[ctg].nswitches;
    }
    INFO(" ");
    INFO("  Total switch errors: %d", switch_errors);
}
