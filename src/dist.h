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

    // supercluster data used to generate this graph
    // query variant data is at sc->callset_vars[QUERY]->{fieldname}[this->qidxs[qni]]
    // truth variant data is at sc->callset_vars[TRUTH]->{fieldname}[this->tidxs[tni]]
    std::shared_ptr<ctgSuperclusters> sc;
	int sc_idx;
    std::string ref;                // for calculating original edit distance
    std::string truth;

    // graph data for each node
    int qnodes;                     // each qvector is of size qnodes
    std::vector<std::string> qseqs; // seq data for each query node (e.g. "ACCCGT")
    std::vector<int> qbegs;         // reference start position
    std::vector<int> qends;         // reference end position
    std::vector<int> qtypes;        // query node TYPE_(REF, INS, SUB, DEL)
    std::vector<int> qidxs;         // store query variant indices (-1 for TYPE_REF)

    int tnodes;                     // each tvector is of size tnodes
    std::vector<std::string> tseqs; // seq data for each truth node (e.g. "ACCCGT")
    std::vector<int> tbegs;         // reference start position
    std::vector<int> ttypes;        // truth node TYPE_(REF, INS, SUB, DEL)
    std::vector<int> tidxs;         // store truth variant indices (-1 for TYPE_REF)

    // set during second pass of graph init
    std::vector< std::vector<int> > qprevs; // directed pointers to prev query nodes
    std::vector< std::vector<int> > qnexts; // directed pointers to next query nodes
    std::vector< std::vector<int> > tprevs; // directed pointers to prev truth nodes
    std::vector< std::vector<int> > tnexts; // directed pointers to next truth nodes

    // constructors
	Graph(std::shared_ptr<ctgSuperclusters> sc, int sc_idx,
			std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hi);

    // methods
    void print();
    int get_truth_pos(int truth_node_idx, int truth_idx);
};

/******************************************************************************/


template <typename T>
bool contains(const std::unordered_set<T> & wave, const T & idx);

template <typename T, typename U>
bool contains(const std::unordered_map<T,U> & wave, const T & idx);

class idx3 { // three indices for alignment cell
public:
    int mi;  // matrix idx
    int qi;  // query idx
    int ri;  // ref idx

    idx3() : mi(0), qi(0), ri(0) {};
    idx3(int q, int r) : mi(0), qi(q), ri(r) {};
    idx3(int m, int q, int r) : mi(m), qi(q), ri(r) {};
    idx3(const idx3 & i2) : mi(i2.mi), qi(i2.qi), ri(i2.ri) {};

    bool operator<(const idx3 & other) const {
        if (this->mi < other.mi) return true;
        else if (this->mi == other.mi && this->qi < other.qi) return true;
        else if (this->mi == other.mi && this->qi == other.qi && this->ri < other.ri) return true;
        return false;
    }
    bool operator==(const idx3 & other) const {
        return this->mi == other.mi && this->qi == other.qi && this->ri == other.ri;
    }
    idx3 & operator=(const idx3 & other) {
        if (this == &other) return *this;
        this->mi = other.mi;
        this->qi = other.qi;
        this->ri = other.ri;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx3> {
        uint64_t operator()(const idx3& x) const noexcept {
            return (uint64_t(x.mi) *73856093 + 0x517cc1b727220a95) ^ 
                   (uint64_t(x.qi) *19349669 + 0xd15f392b3d4704a2) ^ 
                   (uint64_t(x.ri) *83492791);
        }
    };
}

/******************************************************************************/

class idx4 { // four indices for graph cell
public:
    int qni; // query node idx
    int tni; // truth node idx
    int qi;  // query idx
    int ti;  // truth idx

    idx4() : qni(0), tni(0), qi(0), ti(0) {};
    idx4(int qn, int tn, int q, int t) : qni(qn), tni(tn), qi(q), ti(t) {};
    idx4(const idx4 & i2) : qni(i2.qni), tni(i2.tni), qi(i2.qi), ti(i2.ti) {};

    bool operator<(const idx4 & other) const {
        if (this->qni < other.qni) return true;
        if (this->tni < other.tni) return true;
        if (this->qi < other.qi) return true;
        if (this->ti < other.ti) return true;
        return false;
    }
    bool operator==(const idx4 & other) const {
        return this->qni == other.qni && this->tni == other.tni && 
            this->qi == other.qi && this->ti == other.ti;
    }
    bool operator!=(const idx4 & other) const {
        return !(*this == other);
    }
    idx4 & operator=(const idx4 & other) {
        if (this == &other) return *this;
        this->qni = other.qni;
        this->tni = other.tni;
        this->qi = other.qi;
        this->ti = other.ti;
        return *this;
    }
};

namespace std {
    template<> struct hash<idx4> {
        uint64_t operator()(const idx4& x) const noexcept {
            return (uint64_t(x.qni) *73856093 + 0x517cc1b727220a95) ^ 
                   (uint64_t(x.tni) *19349669 + 0xd15f392b3d4704a2) ^ 
                   (uint64_t(x.qi)  *83492791 + 0xc8e13219ab9ab236) ^
                   (uint64_t(x.ti)  *27385201);
        }
    };
}

/******************************************************************************/

std::string generate_str(
        std::shared_ptr<fastaData> ref, 
        std::shared_ptr<ctgVariants> vars, const std::string & ctg,
        int var_idx, int end_idx, int beg_pos, int end_pos, int min_qual=0);

/******************************************************************************/

int calc_ng50(std::vector<int> phase_blocks, size_t total_bases);

/******************************************************************************/

int parse_wfa_path(
        const std::shared_ptr<Graph> graph,
        std::vector< std::vector< std::vector< std::vector<uint32_t> > > > & ptrs
        );


int wfa_calc_prec_recall_aln(
        const std::shared_ptr<Graph> graph,
        std::vector< std::vector< std::vector< std::vector<uint32_t> > > > & ptrs,
        std::vector< std::vector< std::vector< std::vector<int> > > > & offs,
        bool print
        );

int calc_prec_recall_aln(
        const std::shared_ptr<Graph> query_graph,
        std::unordered_map<idx4, idx4> & ptrs,
        bool print = false
        );

void calc_prec_recall(
        const std::shared_ptr<Graph> query_graph,
        std::vector<idx4> & path,
        int truth_hap,
        bool print = false
        );

void evaluate_variants(std::shared_ptr<ctgSuperclusters> sc, int sc_idx,
			std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hi, 
            bool print = false);

/******************************************************************************/

int wf_swg_max_reach(
        const std::string & query, const std::string & truth, 
        std::vector<int> & offs,
        int main_diag, int main_diag_start, int max_score, 
        int x, int o, int e, bool print = false, bool reverse = false
        );

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

/******************************************************************************/

void extract_errors(std::shared_ptr<superclusterData> superclusters,
        std::shared_ptr<variantData> query_ptr_fp,
        std::shared_ptr<variantData> truth_ptr_fn);

/******************************************************************************/

void precision_recall_threads_wrapper(
        std::shared_ptr<superclusterData> clusterdata_ptr,
        std::vector< std::vector< std::vector<int> > > sc_groups);
void precision_recall_wrapper(superclusterData * clusterdata_ptr,
        const std::vector< std::vector< std::vector<int> > > & sc_groups,
        int thread_step, int start, int stop, bool thread2, bool print = false);

int calc_cig_swg_score(
        const std::vector<int> & cigar,
        int sub, int open, int extend);

#endif
