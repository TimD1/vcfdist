#ifndef _DIST_H_
#define _DIST_H_

#include <unordered_set>

#include "fasta.h"
#include "variant.h"
#include "cluster.h"

#define PTR_NONE 0
#define PTR_UP   1
#define PTR_INS  1
#define PTR_LEFT 2
#define PTR_DEL  2
#define PTR_DIAG 4
#define PTR_SUB  8
#define PTR_SWAP 16
#define PTR_LPATH 32
#define PATH      32
#define PTR_RPATH 64
#define MAIN_PATH 96
#define PTR_SYNC 128

#define CALLS1_TRUTH1 0
#define CALLS1_TRUTH2 1
#define CALLS2_TRUTH1 2
#define CALLS2_TRUTH2 3

#define CALLS 0
#define REF 1
#define FALSE 0
#define TRUE 1

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

/******************************************************************************/

class idx {
public:
    int hi;  // hap idx
    int cri; // calls/ref idx
    int ti;  // truth idx

    idx() : hi(0), cri(0), ti(0) {};
    idx(int h, int c, int t) : hi(h), cri(c), ti(t) {};
    idx(const idx & i2) : hi(i2.hi), cri(i2.cri), ti(i2.ti) {};
    bool operator<(const idx & idx2) const {
        if (this->hi < idx2.hi) return true;
        else if (this->hi == idx2.hi && this->cri < idx2.cri) return true;
        else if (this->hi == idx2.hi && this->cri == idx2.cri && this->ti < idx2.ti) return true;
        return false;
    }
    bool operator==(const idx & idx2) const {
        return this->hi == idx2.hi && this->cri == idx2.cri && this->ti == idx2.ti;
    }
    idx & operator=(const idx & other) {
        if (this == &other) return *this;
        this->hi = other.hi;
        this->cri = other.cri;
        this->ti = other.ti;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx> {
        std::uint64_t operator()(const idx& idx1) const noexcept {
            return uint64_t(idx1.hi)<<60 | uint64_t(idx1.cri) << 30 | uint64_t(idx1.ti);
        }
    };
}

bool contains(const std::unordered_set<idx> & wave, const idx & idx);

/******************************************************************************/

void calc_prec_recall_aln(
        std::string calls1, std::string calls2,
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> calls1_ref_ptrs, std::vector<int> ref_calls1_ptrs,
        std::vector<int> calls2_ref_ptrs, std::vector<int> ref_calls2_ptrs,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & ptrs,
        std::vector<int> & pr_calls_ref_end
        );

void get_prec_recall_path_sync(
        std::vector< std::vector<idx> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector<bool> > ref_loc_sync, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector<int> & pr_calls_ref_beg,
        std::vector<int> calls1_ref_ptrs, std::vector<int> ref_calls1_ptrs,
        std::vector<int> calls2_ref_ptrs, std::vector<int> ref_calls2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs
        );

void calc_prec_recall_path(
        std::vector< std::vector<idx> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector<int> calls1_ref_ptrs, std::vector<int> ref_calls1_ptrs,
        std::vector<int> calls2_ref_ptrs, std::vector<int> ref_calls2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs,
        std::vector<int> pr_calls_ref_end
        );

void calc_prec_recall(
        std::shared_ptr<clusterData> clusterdata_ptr, int sc_idx, std::string ctg,
        std::string calls1, std::string calls2, 
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> calls1_ref_ptrs, std::vector<int> ref_calls1_ptrs,
        std::vector<int> calls2_ref_ptrs, std::vector<int> ref_calls2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs,
        std::vector< std::vector<idx> > & path,
        std::vector< std::vector<bool> > & sync,
        std::vector< std::vector<bool> > & edits,
        std::vector< std::vector< std::vector<int> > > & ptrs, 
        std::vector<int> pr_calls_ref_end, int phase, int print
        );

/******************************************************************************/

void alignment_wrapper(std::shared_ptr<clusterData> clusterdata_ptr);

variantData edit_dist_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref);

int calc_vcf_sw_score(
        std::shared_ptr<ctgVariants> vcf, 
        int clust_beg_idx, int clust_end_idx);

void calc_edit_dist_aln(
        std::string calls1, std::string calls2, 
        std::string truth1, std::string truth2,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & offs,
        std::vector< std::vector< std::vector<int> > > & ptrs
        );

#endif
