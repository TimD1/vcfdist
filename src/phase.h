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

    int n = 0;                       // number of phase blocks

    // data, filled out during phase()
    std::vector<int> phase_blocks;   // n+1, variant indices of new phase sets
    std::vector<int> switches;       // switch occurs before this variant
    std::vector<int> flips;          // indices of flipped variants
    int nswitches = 0;               // number of phase switch errors
    int nflips = 0;                  // number of phase flip errors
};

class phaseblockData {
public:
    phaseblockData(std::shared_ptr<superclusterData> clusterdata_ptr); 

    void write_summary_vcf(std::string vcf_fn);
    void write_phasing_summary(int phase_blocks, int switch_errors,
        int flip_errors, int ng50, int s_ngc50, int sf_ngc50);
    void write_switchflips();
    void phase();
    int calculate_ng50(bool break_on_switch = false, bool break_on_flip = false);

    std::shared_ptr<fastaData> ref;
    std::vector<std::string> contigs;
    std::vector<int> lengths;
    std::vector<int> ploidy;
    std::unordered_map<std::string, 
        std::shared_ptr<ctgPhaseblocks> > phase_blocks;
};

#endif
