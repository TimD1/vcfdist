#ifndef _PHASE_H_
#define _PHASE_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "cluster.h"

class ctgPhaseblocks {
public:
    ctgPhaseblocks() {;}

    std::shared_ptr<ctgSuperclusters> ctg_superclusters = nullptr;

    int n = 1;                       // number of phase blocks

    // data, filled out during phase()
    std::vector<int> phase_blocks;   // n+1, indices of superclusters where phasing switches
    std::vector<int> block_states;   // n, phasing for each phase block (keep/swap)
    int nswitches = 0;               // number of phase switch errors
    int nflips = 0;                  // number of supercluster phase flip errors
};

class phaseblockData {
public:
    phaseblockData(std::shared_ptr<superclusterData> clusterdata_ptr); 

    void write_summary_vcf(std::string vcf_fn);
    void phase();
    int calculate_ng50();

    std::shared_ptr<fastaData> ref;
    std::vector<std::string> contigs;
    std::vector<int> lengths;
    std::vector<int> ploidy;
    std::unordered_map<std::string, 
        std::shared_ptr<ctgPhaseblocks> > phase_blocks;
};

#endif
