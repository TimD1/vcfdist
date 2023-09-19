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
    std::vector<int> sc_phasings;    // phasing for each supercluster (keep/swap/either)

    std::shared_ptr<ctgSuperclusters> ctg_superclusters = nullptr;
    int nswitches = 0;               // number of phasing switch errors
    int nerrors = 0;                 // number of phasing errors (whole sc)
    std::vector<int> phase_blocks;   // n+1 indices of superclusters where phasing switches
    std::vector<int> pb_phasings;    // n phasing for each phase block (keep/swap)
};

class phaseData {
public:
    phaseData(std::shared_ptr<superclusterData> clusterdata_ptr); 

    void phase();
    void write_summary_vcf(std::string vcf_fn);

    std::shared_ptr<fastaData> ref;
    std::vector<std::string> contigs;
    std::vector<int> lengths;
    std::vector<int> ploidy;
    std::unordered_map<std::string, 
        std::shared_ptr<ctgPhasings> > ctg_phasings;
};

#endif
