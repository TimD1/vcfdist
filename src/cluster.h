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

    // stores variant info for each callset
    // callset_vars[truth|query]
    std::vector< std::shared_ptr<ctgVariants> > callset_vars;

    // cluster indices of superclusters (not variant indices!)
    // superclusters[truth|query] (size (2, n+1))
    std::vector< std::vector<int> > superclusters;

    // data (length n)
    int n = 0;                                         // number of superclusters
    std::vector<int> begs, ends;                       // ref beg/end pos across haps and truth/query
};

class superclusterData {
public:
    superclusterData(
            std::shared_ptr<variantData> query_ptr,
            std::shared_ptr<variantData> truth_ptr,
            std::shared_ptr<fastaData> ref_ptr);

    void add_callset_vars(int callset, std::vector< std::unordered_map< std::string,
            std::shared_ptr<ctgVariants> > > vars);
    void supercluster();

    // data
    std::vector<std::string> contigs;
    std::vector<std::string> samples;
    std::vector<std::string> filenames;
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
void duplicate_variant_info(std::shared_ptr<superclusterData>);

#endif
