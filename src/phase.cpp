#include <vector>
#include <algorithm>

#include "phase.h"
#include "print.h"
#include "globals.h"

phaseData::phaseData(clusterData* clusters)
{
    // copy contigs
    for (auto ctg : clusters->contigs) {
        this->contigs.push_back(ctg);
        phasings[ctg] = ctgPhasings();
    }

    // for each contig
    for (auto ctg : clusters->contigs) {

        // add pointer to clusters
        this->phasings[ctg].clusters = &(clusters->clusters[ctg]);

        // add phasings
        for(int i = 0; i < clusters->clusters[ctg].n; i++) {
            this->phasings[ctg].phasings.push_back(clusters->clusters[ctg].phase[i]);
            this->phasings[ctg].n++;
        }
    }
    this->phase();
}

void phaseData::phase()
{
    
    // phase each contig separately
    for (auto ctg : this->contigs) {
        std::vector< std::vector<int> > mat(2, std::vector<int>(this->phasings[ctg].n+1));
        std::vector< std::vector<int> > ptrs(2, std::vector<int>(this->phasings[ctg].n));

        // forward pass
        for (int i = 0; i < this->phasings[ctg].n; i++) {

            // determine costs (penalized if this phasing deemed incorrect)
            std::vector<int> costs = {0, 0};
            switch (this->phasings[ctg].phasings[i]) {
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
                    ERROR("Unexpected phase (%d)", this->phasings[ctg].phasings[i]);
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
        int i = this->phasings[ctg].n;
        int phase = 0;
        while (i > 0) {
            if (ptrs[phase][i] == PHASE_PTR_SWAP) {
                phase ^= 1;
                this->phasings[ctg].nswitches++;
                this->phasings[ctg].phase_blocks.push_back(i+1);
            }
            i--;
        }
        this->phasings[ctg].phase_blocks.push_back(0);
        std::reverse(this->phasings[ctg].phase_blocks.begin(),
                this->phasings[ctg].phase_blocks.end());
    }
    

    // print
    for (auto ctg : this->contigs) {

        if (g.print_verbosity >= 2) {
            int s = 0;
            std::vector<std::string> colors {"\033[34m", "\033[32m"};
            int color_idx = 0;
            for (int i = 0; i < this->phasings[ctg].n; i++) {
                if (i == this->phasings[ctg].phase_blocks[s]) {
                    printf("%s", colors[color_idx].data());
                    color_idx ^= 1;
                    s++;
                }
                int ctg_phasing = this->phasings[ctg].phasings[i];
                printf("%s", phase_strs[ctg_phasing].data());
            }
            printf("\033[0m\n");
        }

        if (g.print_verbosity >= 1) {
            printf("phase block cluster indices: ");
            for(int i = 0; i < this->phasings[ctg].nswitches; i++)
                printf("%d ", this->phasings[ctg].phase_blocks[i]);
        }

        INFO("switch errors: %d", this->phasings[ctg].nswitches);
    }

}
