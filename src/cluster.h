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
    ctgSuperclusters();

    // stores variant info for each callset
    // callset_vars[truth|query]
    std::vector< std::shared_ptr<ctgVariants> > callset_vars;

    int get_min_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end);
    int get_max_ref_pos(int qvi_start, int qvi_end, int tvi_start, int tvi_end);
};

class superclusterData {
public:
    superclusterData(
            std::shared_ptr<variantData> query_ptr,
            std::shared_ptr<variantData> truth_ptr,
            std::shared_ptr<fastaData> ref_ptr);

    void add_callset_vars(int callset, std::vector< std::unordered_map< std::string,
            std::shared_ptr<ctgVariants> > > & vars);
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
