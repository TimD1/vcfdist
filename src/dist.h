#ifndef _DIST_H_
#define _DIST_H_

#include <unordered_set>

#include "fasta.h"
#include "variant.h"
#include "cluster.h"

#define PTR_NONE 0  // backtracking pointer flags
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

#define MAT_SUB 0 // three matrices for Smith-Waterman
#define MAT_INS 1
#define MAT_DEL 2
#define MATS    3

#define CALLS1_TRUTH1 0 // four possible alignments of truth and calls haps
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

class idx1 {
public:
    int hi;  // hap idx1
    int cri; // calls/ref idx1
    int ti;  // truth idx1

    idx1() : hi(0), cri(0), ti(0) {};
    idx1(int h, int c, int t) : hi(h), cri(c), ti(t) {};
    idx1(const idx1 & i2) : hi(i2.hi), cri(i2.cri), ti(i2.ti) {};

    bool operator<(const idx1 & other) const {
        if (this->hi < other.hi) return true;
        else if (this->hi == other.hi && this->cri < other.cri) return true;
        else if (this->hi == other.hi && this->cri == other.cri && this->ti < other.ti) return true;
        return false;
    }
    bool operator==(const idx1 & other) const {
        return this->hi == other.hi && this->cri == other.cri && this->ti == other.ti;
    }
    idx1 & operator=(const idx1 & other) {
        if (this == &other) return *this;
        this->hi = other.hi;
        this->cri = other.cri;
        this->ti = other.ti;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx1> {
        std::uint64_t operator()(const idx1& other) const noexcept {
            return uint64_t(other.hi)<<60 | uint64_t(other.cri) << 30 | uint64_t(other.ti);
        }
    };
}

template <typename T>
bool contains(const std::unordered_set<T> & wave, const T & idx);

class idx2 {
public:
    int mi;  // matrix idx2
    int ci;  // calls idx2
    int ri;  // ref idx2

    idx2() : mi(0), ci(0), ri(0) {};
    idx2(int c, int r) : mi(0), ci(c), ri(r) {};
    idx2(int m, int c, int r) : mi(m), ci(c), ri(r) {};
    idx2(const idx2 & i2) : mi(i2.mi), ci(i2.ci), ri(i2.ri) {};

    bool operator<(const idx2 & other) const {
        if (this->mi < other.mi) return true;
        else if (this->mi == other.mi && this->ci < other.ci) return true;
        else if (this->mi == other.mi && this->ci == other.ci && this->ri < other.ri) return true;
        return false;
    }
    bool operator==(const idx2 & other) const {
        return this->mi == other.mi && this->ci == other.ci && this->ri == other.ri;
    }
    idx2 & operator=(const idx2 & other) {
        if (this == &other) return *this;
        this->mi = other.mi;
        this->ci = other.ci;
        this->ri = other.ri;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx2> {
        std::uint64_t operator()(const idx2& other) const noexcept {
            return uint64_t(other.mi)<<60 | uint64_t(other.ci) << 30 | uint64_t(other.ri);
        }
    };
}

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
        std::vector< std::vector<idx1> > & path, 
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
        std::vector< std::vector<idx1> > & path, 
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
        std::vector< std::vector<idx1> > & path,
        std::vector< std::vector<bool> > & sync,
        std::vector< std::vector<bool> > & edits,
        std::vector< std::vector< std::vector<int> > > & ptrs, 
        std::vector<int> pr_calls_ref_end, int phase, int print
        );

/******************************************************************************/

int sw_max_reach_ref(std::string calls, std::string ref, 
        std::vector<int> calls_ref_ptrs, std::vector<int> ref_calls_ptrs,
        int score);

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
