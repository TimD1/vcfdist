#ifndef _DIST_H_
#define _DIST_H_

#include "fasta.h"
#include "variant.h"
#include "cluster.h"

#define PTR_NONE 0
#define PTR_UP   1
#define PTR_LEFT 2
#define PTR_BOTH 3
#define PTR_DIAG 4
#define PTR_SWAP 8
#define PTR_LPATH 16
#define PTR_RPATH 32
#define PTR_SUB 64
#define PTR_DONE 128
#define PTR_NEXT 256

#define CALLS1_TRUTH1 0
#define CALLS1_TRUTH2 1
#define CALLS2_TRUTH1 2
#define CALLS2_TRUTH2 3

#define CALLS 0
#define REF 1

#define ERRTYPE_TP 0 // true positive
#define ERRTYPE_FP 1 // false positive
#define ERRTYPE_FN 2 // false negative
#define ERRTYPE_PP 3 // partial positive (reduces ED, but not TP)
#define ERRTYPE_PE 4 // phase error (0|1 -> 1|0)
#define ERRTYPE_GE 5 // genotype error (0|1 -> 1|1)
#define ERRTYPE_UN 6 // unknown
#define ERRTYPES 7

void generate_ptrs_strs(
        std::string & hap1, std::string & hap2,          // actual strings 
        std::string & hap1_str, std::string & hap2_str,  // colored debug strs
        std::vector<int> & hap1_ptrs, std::vector<int> & hap2_ptrs,
        std::shared_ptr<ctgVariants> hap1_vars, std::shared_ptr<ctgVariants> hap2_vars,
        size_t hap1_clust_beg_idx, size_t hap2_clust_beg_idx,
        size_t hap1_clust_end_idx, size_t hap2_clust_end_idx,
        int beg_pos, int end_pos, std::shared_ptr<fastaData> ref, std::string ctg
        );

variantData edit_dist_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref_ptr);

void alignment_wrapper(std::shared_ptr<clusterData> clusterdata_ptr);

void edit_dist(
        std::string calls1, std::string calls2, 
        std::vector<int> calls1_ptrs,
        std::vector<int> calls2_ptrs, 
        std::string truth1, std::string truth2,
        std::vector<int> truth1_ptrs,
        std::vector<int> truth2_ptrs, 
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & offs,
        std::vector< std::vector< std::vector<int> > > & ptrs
        );

#endif
