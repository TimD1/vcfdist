#ifndef _DIST_H_
#define _DIST_H_

#include <unordered_set>
#include <unordered_map>

#include "fasta.h"
#include "variant.h"
#include "cluster.h"
#include "defs.h"
#include "edit.h"

void generate_ptrs_strs(
        std::string & query_str, std::string & ref_str,
        std::vector< std::vector<int> > & query_ptrs, 
        std::vector< std::vector<int> > & ref_ptrs,
        std::shared_ptr<ctgVariants> query_vars,
        size_t query_clust_beg_idx, size_t ref_clust_beg_idx,
        int beg_pos, int end_pos, std::shared_ptr<fastaData> ref, 
        const std::string & ctg
        );

void reverse_ptrs_strs(
        std::string & query, std::string & ref,
        std::vector< std::vector<int> > & query_ptrs, 
        std::vector< std::vector<int> > & ref_ptrs
        );

/******************************************************************************/

class idx1 {
public:
    int hi;  // hap idx1
    int qri; // query/ref idx1
    int ti;  // truth idx1

    idx1() : hi(0), qri(0), ti(0) {};
    idx1(int h, int q, int t) : hi(h), qri(q), ti(t) {};
    idx1(const idx1 & i2) : hi(i2.hi), qri(i2.qri), ti(i2.ti) {};

    bool operator<(const idx1 & other) const {
        if (this->hi < other.hi) return true;
        else if (this->hi == other.hi && this->qri < other.qri) return true;
        else if (this->hi == other.hi && this->qri == other.qri && this->ti < other.ti) return true;
        return false;
    }
    bool operator==(const idx1 & other) const {
        return this->hi == other.hi && this->qri == other.qri && this->ti == other.ti;
    }
    idx1 & operator=(const idx1 & other) {
        if (this == &other) return *this;
        this->hi = other.hi;
        this->qri = other.qri;
        this->ti = other.ti;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx1> {
        std::uint64_t operator()(const idx1& other) const noexcept {
            return uint64_t(other.hi)<<60 | uint64_t(other.qri) << 30 | uint64_t(other.ti);
        }
    };
}

template <typename T>
bool contains(const std::unordered_set<T> & wave, const T & idx);

class idx2 {
public:
    int mi;  // matrix idx2
    int qi;  // query idx2
    int ri;  // ref idx2

    idx2() : mi(0), qi(0), ri(0) {};
    idx2(int q, int r) : mi(0), qi(q), ri(r) {};
    idx2(int m, int q, int r) : mi(m), qi(q), ri(r) {};
    idx2(const idx2 & i2) : mi(i2.mi), qi(i2.qi), ri(i2.ri) {};

    bool operator<(const idx2 & other) const {
        if (this->mi < other.mi) return true;
        else if (this->mi == other.mi && this->qi < other.qi) return true;
        else if (this->mi == other.mi && this->qi == other.qi && this->ri < other.ri) return true;
        return false;
    }
    bool operator==(const idx2 & other) const {
        return this->mi == other.mi && this->qi == other.qi && this->ri == other.ri;
    }
    idx2 & operator=(const idx2 & other) {
        if (this == &other) return *this;
        this->mi = other.mi;
        this->qi = other.qi;
        this->ri = other.ri;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx2> {
        std::uint64_t operator()(const idx2& other) const noexcept {
            return uint64_t(other.mi)<<60 | uint64_t(other.qi) << 30 | uint64_t(other.ri);
        }
    };
}

/******************************************************************************/

void calc_prec_recall_aln(
        const std::string & query1, const std::string & query2,
        const std::string & truth1, const std::string & truth2, 
        const std::string & ref,
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs,
        std::vector< std::unordered_map<idx1, idx1> > & pred_map,
        std::vector<int> & pr_query_ref_end, bool print
        );

void get_prec_recall_path_sync(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        const std::vector< std::vector<bool> > & ref_loc_sync, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_beg, bool print
        );

void calc_prec_recall_path(
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        std::vector< std::vector< std::vector<int> > > & path_scores, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector< std::unordered_map<idx1, idx1> > & pred_map,
        const std::vector<int> & pr_query_ref_end, bool print
        );

void calc_prec_recall(
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, 
        const std::string & ctg, const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        const std::vector< std::vector<idx1> > & path,
        const std::vector< std::vector<bool> > & sync,
        const std::vector< std::vector<bool> > & edits,
        const std::vector< std::vector< std::vector<int> > > & aln_ptrs, 
        const std::vector< std::vector< std::vector<int> > > & path_ptrs, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_end, int phase, bool print
        );

/******************************************************************************/

int sw_max_reach(
        const std::string & query, 
        const std::string & ref, 
        const std::vector< std::vector<int> > & query_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query_ptrs,
        int sub, int open, int extend,
        int score, bool print,
        bool reverse = false,
        int ref_section = -1);

std::unique_ptr<variantData> sw_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref,
        int sub, int open, int extend, int callset, bool print = false);

std::unordered_map<idx2,idx2> sw_align(
        const std::string & query, 
        const std::string & truth,
        int sub, int open, int extend, bool print);

std::vector<int> sw_backtrack(
        const std::string & query,
        const std::string & truth,
        const std::unordered_map<idx2, idx2> & ptrs, bool print);

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
int calc_cig_sw_score(
        const std::vector<int> & cigar,
        int sub, int open, int extend);

void calc_edit_dist_aln(
        const std::string & query1, 
        const std::string & query2, 
        const std::string & truth1, 
        const std::string & truth2,
        std::vector<int> & s, 
        std::vector< std::vector< std::vector<int> > > & offs,
        std::vector< std::vector< std::vector<int> > > & ptrs
        );

#endif
