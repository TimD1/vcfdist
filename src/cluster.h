#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "fasta.h"
#include "variant.h"
#include "defs.h"

class ctgSuperclusters {
public:

    // constructor
    ctgSuperclusters();

    // add supercluster info
    void add_supercluster(
           std::vector<int> brks,
           int beg, int end);
    void set_phase(
           int sc_idx, int phase, 
           int orig_phase_dist, 
           int swap_phase_dist);

    // stores variant info for each contig
    // ctg_variants[truth/query][hap]
    std::vector< std::vector< std::shared_ptr<ctgVariants> > > ctg_variants;

    // cluster indices of superclusters (not variant indices!)
    // superclusters[truth/query][hap] = n+1
    std::vector< std::vector< std::vector<int> > > superclusters;

    // data (length n)
    int n = 0;                                         // number of superclusters
    std::vector<int> begs, ends;                       // ref beg/end pos across haps and truth/query
    std::vector<int> sc_phase;                         // keep/swap/unknown, from alignment
    std::vector<int> pb_phase;                         // keep/swap, from phasing algorithm
    std::vector<int> phase_sets;                       // input, from first variant in sc
    std::vector<int> orig_phase_dist, swap_phase_dist; // alignment distances
};

class superclusterData {
public:
    superclusterData(
            std::shared_ptr<variantData> query_ptr,
            std::shared_ptr<variantData> truth_ptr,
            std::shared_ptr<fastaData> ref_ptr);

    void supercluster();
    void transfer_phase_sets();

    // data
    std::vector<std::string> contigs;
    std::vector<int> lengths;
    std::vector<int> ploidy;
    std::unordered_map<std::string, 
        std::shared_ptr<ctgSuperclusters> > superclusters;
    std::shared_ptr<fastaData> ref;
};


// for single haplotype clustering (one VCF)
void simple_cluster(std::shared_ptr<variantData> vcf, int callset);
void wf_swg_cluster(variantData * vcf, int ctg_idx, int hap,
        int sub, int open, int extend);
std::vector< std::vector< std::vector<int> > > 
        sort_superclusters(std::shared_ptr<superclusterData>);

#endif
