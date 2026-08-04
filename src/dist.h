/**
 * @file dist.h
 * @brief Graph-based alignment and precision/recall evaluation declarations.
 */
#ifndef _DIST_H_
#define _DIST_H_

#include <unordered_set>
#include <unordered_map>

#include "fasta.h"
#include "variant.h"
#include "cluster.h"
#include "defs.h"

/**
 * @class Graph
 * @brief Store graph information for truth and query variants.
 */
class Graph {
public:

    /** @brief Store supercluster data used to generate this graph
        query variant data is at sc->callset_vars[QUERY]->{fieldname}[this->qidxs[qni]]
        truth variant data is at sc->callset_vars[TRUTH]->{fieldname}[this->tidxs[tni]]
    */
    std::shared_ptr<ctgSuperclusters> sc;
	int sc_idx;                     ///< supercluster index for this graph
    std::string ref;                ///< reference sequence, for calculating original edit distance
    std::string truth;              ///< truth sequence

    // graph data for each node
    int qnodes;                     ///< each qvector is of size qnodes
    std::vector<std::string> qseqs; ///< seq data for each query node (e.g. "ACCCGT")
    std::vector<int> qbegs;         ///< reference start position
    std::vector<int> qends;         ///< reference end position
    std::vector<int> qtypes;        ///< query node TYPE_(REF, INS, SUB, DEL)
    std::vector<int> qidxs;         ///< store query variant indices (-1 for TYPE_REF)

    int tnodes;                     ///< each tvector is of size tnodes
    std::vector<std::string> tseqs; ///< seq data for each truth node (e.g. "ACCCGT")
    std::vector<int> tbegs;         ///< reference start position
    std::vector<int> ttypes;        ///< truth node TYPE_(REF, INS, SUB, DEL)
    std::vector<int> tidxs;         ///< store truth variant indices (-1 for TYPE_REF)

    // set during second pass of graph init
    std::vector< std::vector<int> > qprevs; ///< directed pointers to prev query nodes
    std::vector< std::vector<int> > qnexts; ///< directed pointers to next query nodes
    std::vector< std::vector<int> > tprevs; ///< directed pointers to prev truth nodes
    std::vector< std::vector<int> > tnexts; ///< directed pointers to next truth nodes

    /** @brief Constructs alignment graph from supercluster variants and reference sequence. */
    Graph(std::shared_ptr<ctgSuperclusters> sc, int sc_idx,
            std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hi);

    /** @brief Prints graph node sequences and connectivity to console for debugging. */
    void print();

    /** @brief Maps a truth graph node and variant index to its 0-based reference position. */
    int get_truth_pos(int truth_node_idx, int truth_idx);
};

/**************************************************************************************************/

/**
 * @class idx4
 * @brief Store the indices for a graph cell.
 *
 * Node indices define the current submatrix, remaining indices define the current cell.
 */
class idx4 { // four indices for graph cell
public:
    int qni; ///< query node idx
    int tni; ///< truth node idx
    int qi;  ///< query idx
    int ti;  ///< truth idx

    /** @brief Default constructor; initializes all indices to zero. */
    idx4() : qni(0), tni(0), qi(0), ti(0) {};

    /**
     * @brief Constructs idx4 with specified node and position indices.
     * @param[in] qn Query node index
     * @param[in] tn Truth node index
     * @param[in] q Query position index within node
     * @param[in] t Truth position index within node
     */
    idx4(int qn, int tn, int q, int t) : qni(qn), tni(tn), qi(q), ti(t) {};

    /**
     * @brief Copy constructor.
     * @param[in] i2 Source idx4 to copy
     */
    idx4(const idx4 & i2) : qni(i2.qni), tni(i2.tni), qi(i2.qi), ti(i2.ti) {};

    /**
     * @brief Lexicographic less-than comparison for use in ordered containers.
     * @return True if this cell precedes other in (qni, tni, qi, ti) order
     */
    bool operator<(const idx4 & other) const {
        if (this->qni != other.qni) return this->qni < other.qni;
        if (this->tni != other.tni) return this->tni < other.tni;
        if (this->qi  != other.qi)  return this->qi  < other.qi;
        return this->ti < other.ti;
    }

    /**
     * @brief Equality comparison; true when all four indices match.
     * @return True if all fields are equal
     */
    bool operator==(const idx4 & other) const {
        return this->qni == other.qni && this->tni == other.tni &&
            this->qi == other.qi && this->ti == other.ti;
    }

    /**
     * @brief Inequality comparison.
     * @return True if any field differs
     */
    bool operator!=(const idx4 & other) const {
        return !(*this == other);
    }

    /**
     * @brief Copy assignment operator.
     * @return Reference to this after assignment
     */
    idx4 & operator=(const idx4 & other) {
        if (this == &other) return *this;
        this->qni = other.qni;
        this->tni = other.tni;
        this->qi = other.qi;
        this->ti = other.ti;
        return *this;
    }
};

/**
 * @brief Hash specialization enabling idx4 as unordered_map key.
 */
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

/**************************************************************************************************/

/** @brief Returns true if a hash set contains the given value. */
template <typename T>
bool contains(const std::unordered_set<T> & hash_set, const T & value) {
    return hash_set.find(value) != hash_set.end();
}

/** @brief Generates a haplotype sequence string by applying variants to the reference. */
std::string generate_str(
        std::shared_ptr<fastaData> ref,
        std::shared_ptr<ctgVariants> vars, const std::string & ctg,
        int var_idx, int end_idx, int beg_pos, int end_pos, int min_qual=0);

/**************************************************************************************************/

/** @brief Calculates the NG50 statistic for a set of phase block lengths. */
int calc_ng50(std::vector<int> phase_blocks, size_t total_bases);

/**************************************************************************************************/

/** @brief Evaluates query variants against truth for one supercluster and haplotype combination. */
void evaluate_variants(std::shared_ptr<ctgSuperclusters> sc, int sc_idx,
			std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hi,
            bool print = false);

/** @brief Runs graph-based alignment and returns the optimal alignment score. */
int calc_prec_recall_aln(
        const std::shared_ptr<Graph> query_graph,
        std::unordered_map<idx4, idx4> & ptrs,
        bool print = false
        );

/** @brief Assigns TP/FP/FN error types and credit scores to variants from an alignment. */
void calc_prec_recall(
        const std::shared_ptr<Graph> query_graph,
        const std::unordered_map<idx4, idx4> & ptrs,
        int truth_hap,
        bool print = false
        );

/** @brief Launches threaded precision/recall evaluation across all superclusters. */
void precision_recall_threads_wrapper(
        std::shared_ptr<superclusterData> clusterdata_ptr,
        std::vector< std::vector< std::vector<int> > > sc_groups);

/** @brief Evaluates a subset of superclusters within a single thread. */
void precision_recall_wrapper(superclusterData * clusterdata_ptr,
        const std::vector< std::vector< std::vector<int> > > & sc_groups,
        int thread_step, int start, int stop, bool thread2, bool print = false);

/**************************************************************************************************/

/** @brief Extends a wavefront diagonal to its maximum reach using Smith-Waterman gap scoring. */
int wf_swg_max_reach(
        const std::string & query, const std::string & truth,
        std::vector<int> & offs,
        int main_diag, int main_diag_start, int max_score,
        int x, int o, int e, bool print = false, bool reverse = false
        );

/** @brief Computes Smith-Waterman gap-affine alignment score between two sequences. */
void wf_swg_align(
        const std::string & query,
        const std::string & truth,
        int & score, int sub, int open, int extend, bool print = false);

/** @brief Computes the edit distance between two sequences using wavefront alignment. */
void wf_ed(const std::string & query, const std::string & truth, int & score, bool print = false);

#endif
