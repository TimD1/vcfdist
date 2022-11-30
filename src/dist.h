#ifndef _DIST_H_
#define _DIST_H_

#include <unordered_set>

#include "fasta.h"
#include "variant.h"
#include "cluster.h"
#include "defs.h"
#include "edit.h"

void generate_ptrs_strs(
        std::string & hap1, std::string & hap2,
        std::vector<int> & hap1_ptrs, std::vector<int> & hap2_ptrs,
        std::shared_ptr<ctgVariants> hap1_vars, std::shared_ptr<ctgVariants> hap2_vars,
        size_t hap1_clust_beg_idx, size_t hap2_clust_beg_idx,
        size_t hap1_clust_end_idx, size_t hap2_clust_end_idx,
        int beg_pos, int end_pos, std::shared_ptr<fastaData> ref, std::string ctg
        );

void reverse_ptrs_strs(std::string & query, std::string & ref,
        std::vector<int> & query_ptrs, std::vector<int> & ref_ptrs);

/******************************************************************************/

class idx1 {
public:
    int hi;  // hap idx1
    int cri; // query/ref idx1
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
    int ci;  // query idx2
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
        std::string query1, std::string query2,
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & ptrs,
        std::vector<int> & pr_query_ref_end
        );

void get_prec_recall_path_sync(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector<bool> > ref_loc_sync, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector<int> & pr_query_ref_beg,
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs
        );

void calc_prec_recall_path(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs,
        std::vector<int> pr_query_ref_end
        );

void calc_prec_recall(
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, std::string ctg,
        std::string query1, std::string query2, 
        std::string truth1, std::string truth2, std::string ref,
        std::vector<int> query1_ref_ptrs, std::vector<int> ref_query1_ptrs,
        std::vector<int> query2_ref_ptrs, std::vector<int> ref_query2_ptrs,
        std::vector<int> truth1_ref_ptrs, std::vector<int> truth2_ref_ptrs,
        std::vector< std::vector<idx1> > & path,
        std::vector< std::vector<bool> > & sync,
        std::vector< std::vector<bool> > & edits,
        std::vector< std::vector< std::vector<int> > > & ptrs, 
        std::vector<int> pr_query_ref_end, int phase, int print
        );

/******************************************************************************/

int sw_max_reach(
        std::string query, 
        std::string ref, 
        std::vector<int> query_ref_ptrs, 
        std::vector<int> ref_query_ptrs,
        int sub, int open, int extend,
        int score, 
        bool reverse=false);

std::unique_ptr<variantData> sw_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref,
        int sub, int open, int extend);

std::unordered_map<idx2,idx2> sw_align(
        const std::string & query, 
        const std::string & truth,
        int sub, int open, int extend);

std::vector<int> sw_backtrack(
        const std::string & query,
        const std::string & truth,
        const std::unordered_map<idx2, idx2> & ptrs);

int count_dist(const std::vector<int> & cigar);

/******************************************************************************/

editData alignment_wrapper(
        std::shared_ptr<superclusterData> clusterdata_ptr);

variantData edit_dist_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref);

int calc_vcf_sw_score(
        std::shared_ptr<ctgVariants> vcf, 
        int clust_beg_idx, 
        int clust_end_idx,
        int sub, int open, int extend);

void calc_edit_dist_aln(
        std::string query1, 
        std::string query2, 
        std::string truth1, 
        std::string truth2,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & offs,
        std::vector< std::vector< std::vector<int> > > & ptrs
        );

#endif
