#ifndef _PHASE_H_
#define _PHASE_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "cluster.h"

class ctgPhasings {
public:
    ctgPhasings() {;}

    int n = 0;                       // number of superclusters
    std::vector<int> phasings;       // phasing for each supercluster (keep/swap)

    std::shared_ptr<ctgSuperclusters> ctg_superclusters = nullptr;
    int nswitches = 0;               // number of phasing errors
    std::vector<int> phase_blocks;   // indices of superclusters where phasing switches
};

class phaseData {
public:
    phaseData(std::shared_ptr<superclusterData> clusterdata_ptr); 

    void phase();

    std::vector<std::string> contigs;
    std::unordered_map<std::string, ctgPhasings> ctg_phasings;
};

#endif
