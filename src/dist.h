#ifndef _DIST_H_
#define _DIST_H_

#include <unordered_set>
#include <unordered_map>

#include "fasta.h"
#include "variant.h"
#include "cluster.h"
#include "defs.h"
#include "edit.h"


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
        std::uint64_t operator()(const idx1& x) const noexcept {
            return (uint64_t(x.hi) *73856093 + 0x517cc1b727220a95) ^ 
                   (uint64_t(x.qri)*19349669 + 0xd15f392b3d4704a2) ^ 
                   (uint64_t(x.ti) *83492791);
        }
    };
}

template <typename T>
bool contains(const std::unordered_set<T> & wave, const T & idx);

template <typename T, typename U>
bool contains(const std::unordered_map<T,U> & wave, const T & idx);

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
        std::uint64_t operator()(const idx2& x) const noexcept {
            return (uint64_t(x.mi) *73856093 + 0x517cc1b727220a95) ^ 
                   (uint64_t(x.qi) *19349669 + 0xd15f392b3d4704a2) ^ 
                   (uint64_t(x.ri) *83492791);
        }
    };
}

/******************************************************************************/

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

std::string generate_str(
        std::shared_ptr<fastaData> ref, 
        std::shared_ptr<ctgVariants> vars, const std::string & ctg,
        int var_idx, int end_idx, int beg_pos, int end_pos, int min_qual=0);

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
        std::vector< std::vector< std::vector<uint8_t> > > & ptrs, 
        std::unordered_map<idx1, idx1> & swap_pred_map,
        std::vector<int> & pr_query_ref_end, bool print
        );

void calc_prec_recall_path(
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<uint8_t> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<uint8_t> > > & path_ptrs, 
        std::vector< std::vector< std::vector<int16_t> > > & path_scores, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::unordered_map<idx1,idx1> & swap_pred_map,
        const std::vector<int> & pr_query_ref_end, int phase, bool print
        );

void get_prec_recall_path_sync(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<uint8_t> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<uint8_t> > > & path_ptrs, 
        const std::vector< std::vector<bool> > & ref_loc_sync, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_beg, int phase, bool print
        );

void calc_prec_recall(
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, 
        const std::string & ctg, 
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        const std::vector< std::vector<idx1> > & path,
        const std::vector< std::vector<bool> > & sync,
        const std::vector< std::vector<bool> > & edits,
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_end, int phase, 
        bool print
        );


/******************************************************************************/

int swg_max_reach(
        const std::string & query, 
        const std::string & ref, 
        const std::vector< std::vector<int> > & query_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query_ptrs,
        int sub, int open, int extend,
        int score, bool print,
        bool reverse = false,
        int ref_section = -1);

int wf_swg_max_reach(
        const std::string & query, const std::string & truth, 
        std::vector<int> & offs,
        int main_diag, int main_diag_start, int max_score, 
        int x, int o, int e, bool print = false, bool reverse = false
        );

std::unique_ptr<variantData> wf_swg_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref,
        int sub, int open, int extend, int callset, bool print = false);

void wf_swg_align(
        const std::string & query, 
        const std::string & truth,
        const std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        const std::vector< std::vector< std::vector<int> > > & offs,
        int & s, int sub, int open, int extend, bool print = false);

std::vector<int> wf_swg_backtrack(
        const std::string & query, 
        const std::string & truth,
        const std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        const std::vector< std::vector< std::vector<int> > > & offs,
        int s, int sub, int open, int extend, bool print = false);

void wf_ed(
        const std::string & query, const std::string & truth, int & s, 
        std::vector< std::vector<int> > & offs,
        std::vector< std::vector<int> > & ptrs);


int count_dist(const std::vector<int> & cigar);

/******************************************************************************/

editData alignment_wrapper(
        std::shared_ptr<superclusterData> clusterdata_ptr);

int calc_vcf_swg_score(
        std::shared_ptr<ctgVariants> vcf, 
        int clust_beg_idx, 
        int clust_end_idx,
        int sub, int open, int extend);

int calc_cig_swg_score(
        const std::vector<int> & cigar,
        int sub, int open, int extend);

#endif
