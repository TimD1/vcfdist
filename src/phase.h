#ifndef _PHASE_H_
#define _PHASE_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "cluster.h"

class ctgPhasings {
public:

    // constructor
    ctgPhasings() {;}

    int n = 0;                       // number of superclusters
    std::vector<int> phasings;       // phasing for each supercluster (keep/swap)

    ctgClusters* clusters = nullptr;
    int nswitches = 0;               // number of phasing errors
    std::vector<int> phase_blocks;   // indices of superclusters where phasing switches
};

class phaseData {
public:
    //constructor
    phaseData(clusterData* clusters); 

    // member
    void phase();

    // data
    std::vector<std::string> contigs;
    std::unordered_map<std::string, ctgPhasings> phasings;
};

#endif
