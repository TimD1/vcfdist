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

    int n = 0;
    std::vector<int> phasings;
    ctgClusters* clusters = nullptr;
    int nswitches = 0;
    std::vector<int> phase_blocks;
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
