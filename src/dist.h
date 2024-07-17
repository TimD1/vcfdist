#ifndef _DIST_H_
#define _DIST_H_

#include <unordered_set>
#include <unordered_map>

#include "fasta.h"
#include "variant.h"
#include "cluster.h"
#include "defs.h"

class Graph {
public:

    // information about this graph (generated for a supercluster)
    std::shared_ptr<ctgSuperclusters> sc;
	int sc_idx;

    // data for each node (all size n)
    int n;
    // NOTE: some of this information is unnecessarily duplicated, since most can be accessed 
    // through idxs sc->callset_vars[QUERY]->{fieldname}[this->idxs[i]]
    std::vector<std::string> seqs;         // seq data for each node (e.g. "ACCCGT")
    std::vector<int> begs;                 // node reference start positions (relative to sc->begs[sc_idx])
    std::vector<int> ends;                 // node reference end positions (relative to sc->begs[sc_idx])
    std::vector<int> types;                // node TYPE_(REF, INS, SUB, DEL)
    // for TYPE_REF, idxs is -1
    std::vector<int> idxs;                 // store variant indices for checking per-variant data
    // set during second pass of graph init
    std::vector< std::vector<int> > prevs; // directed pointers to previous nodes
    std::vector< std::vector<int> > nexts; // directed pointers to next nodes

    // constructors
	Graph(std::shared_ptr<ctgSuperclusters> sc, int sc_idx,
			std::shared_ptr<fastaData> ref, std::string ctg);

    // methods
    void print();
};

/******************************************************************************/


template <typename T>
bool contains(const std::unordered_set<T> & wave, const T & idx);

template <typename T, typename U>
bool contains(const std::unordered_map<T,U> & wave, const T & idx);

class idx {
public:
    int mi;  // matrix idx
    int qi;  // query idx
    int ri;  // ref idx

    idx() : mi(0), qi(0), ri(0) {};
    idx(int q, int r) : mi(0), qi(q), ri(r) {};
    idx(int m, int q, int r) : mi(m), qi(q), ri(r) {};
    idx(const idx & i2) : mi(i2.mi), qi(i2.qi), ri(i2.ri) {};

    bool operator<(const idx & other) const {
        if (this->mi < other.mi) return true;
        else if (this->mi == other.mi && this->qi < other.qi) return true;
        else if (this->mi == other.mi && this->qi == other.qi && this->ri < other.ri) return true;
        return false;
    }
    bool operator==(const idx & other) const {
        return this->mi == other.mi && this->qi == other.qi && this->ri == other.ri;
    }
    idx & operator=(const idx & other) {
        if (this == &other) return *this;
        this->mi = other.mi;
        this->qi = other.qi;
        this->ri = other.ri;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx> {
        uint64_t operator()(const idx& x) const noexcept {
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
        const std::string & ctg, int hap
        );

std::string generate_str(
        std::shared_ptr<fastaData> ref, 
        std::shared_ptr<ctgVariants> vars, const std::string & ctg,
        int var_idx, int end_idx, int beg_pos, int end_pos, int min_qual=0);

/******************************************************************************/

int calc_ng50(std::vector<int> phase_blocks, size_t total_bases);

/******************************************************************************/

int calc_prec_recall_aln(
        const std::shared_ptr<Graph> query_graph,
        const std::string & truth,
        const std::string & ref,
        std::unordered_map<idx, idx> & ptrs,
        bool print
        );

void calc_prec_recall_path(
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        std::vector< std::vector<idx> > & path, 
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
        const std::vector< std::shared_ptr< std::unordered_map<idx,idx> > > & swap_pred_maps,
        const std::vector<int> & pr_query_ref_end, bool print
        );

void get_prec_recall_path_sync(
        std::vector< std::vector<idx> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::vector< std::vector< std::vector<uint8_t> > > & aln_ptrs, 
        std::vector< std::vector< std::vector<uint8_t> > > & path_ptrs, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_beg, bool print
        );

void calc_prec_recall(
        superclusterData * clusterdata_ptr, int sc_idx, 
        const std::string & ctg, 
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        const std::vector< std::vector<idx> > & path,
        const std::vector< std::vector<bool> > & sync,
        const std::vector< std::vector<bool> > & edits,
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_end, bool print
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

std::shared_ptr<variantData> wf_swg_realign(
        std::shared_ptr<variantData> vcf, 
        std::shared_ptr<fastaData> ref,
        int sub, int open, int extend, int callset, bool print = false);

void wf_swg_align(
        const std::string & query, 
        const std::string & truth,
        std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        std::vector< std::vector< std::vector<int> > > & offs,
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
        std::vector< std::vector<int> > & ptrs, bool print = false);


int count_dist(const std::vector<int> & cigar);

/******************************************************************************/

void precision_recall_threads_wrapper(
        std::shared_ptr<superclusterData> clusterdata_ptr,
        std::vector< std::vector< std::vector<int> > > sc_groups);
void precision_recall_wrapper(superclusterData * clusterdata_ptr,
        const std::vector< std::vector< std::vector<int> > > & sc_groups,
        int thread_step, int start, int stop, bool thread4);

int calc_vcf_swg_score(
        std::shared_ptr<ctgVariants> vcf, 
        int clust_beg_idx, 
        int clust_end_idx,
        int sub, int open, int extend);

int calc_cig_swg_score(
        const std::vector<int> & cigar,
        int sub, int open, int extend);

#endif
