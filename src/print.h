/**
 * @file print.h
 * @brief Terminal color utilities, quality scoring, and output writing functions.
 */
#ifndef _PRINT_H_
#define _PRINT_H_

#include <memory>
#include <string>
#include <vector>

#include "globals.h"
#include "phase.h"
#include "defs.h"
#include "dist.h"

/** @brief Wraps an integer in ANSI green codes for terminal output. */
std::string GREEN(int i);

/** @brief Wraps a char in ANSI green codes for terminal output. */
std::string GREEN(char c);

/** @brief Wraps a string in ANSI green codes for terminal output. */
std::string GREEN(const std::string & str);

/** @brief Wraps an integer in ANSI red codes for terminal output. */
std::string RED(int i);

/** @brief Wraps a char in ANSI red codes for terminal output. */
std::string RED(char c);

/** @brief Wraps a string in ANSI red codes for terminal output. */
std::string RED(const std::string & str);

/** @brief Wraps an integer in ANSI blue codes for terminal output. */
std::string BLUE(int i);

/** @brief Wraps a char in ANSI blue codes for terminal output. */
std::string BLUE(char c);

/** @brief Wraps a string in ANSI blue codes for terminal output. */
std::string BLUE(const std::string & str);

/** @brief Wraps an integer in ANSI yellow codes for terminal output. */
std::string YELLOW(int i);

/** @brief Wraps a char in ANSI yellow codes for terminal output. */
std::string YELLOW(char c);

/** @brief Wraps a string in ANSI yellow codes for terminal output. */
std::string YELLOW(const std::string & str);

/** @brief Wraps an integer in ANSI purple codes for terminal output. */
std::string PURPLE(int i);

/** @brief Wraps a char in ANSI purple codes for terminal output. */
std::string PURPLE(char c);

/** @brief Wraps a string in ANSI purple codes for terminal output. */
std::string PURPLE(const std::string & str);

/** @brief Converts error probability to Phred quality score. */
float qscore(double p_error);

/** @brief Returns string representation of a graph cell's alignment pointer. */
std::string get_ptr_repr(idx4 cell, const std::unordered_map<idx4,idx4> & ptrs);

/** @brief Prints alignment matrix for graph-based alignment (debug utility). */
void print_graph_ptrs(const std::shared_ptr<Graph> graph,
        const std::unordered_map<idx4,idx4> & ptrs);

/** @brief Prints WFA substitution/insertion/deletion matrices for debugging. */
void print_wfa_ptrs(
        const std::string & query,
        const std::string & truth,
        int s,
        const std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        const std::vector< std::vector< std::vector<int> > > & offs);

/**
 * @struct pr_counts
 * @brief Variant counts for each callset, error type, and quality threshold.
 */
struct pr_counts {
    /// query counts; ax0: SNP/INDEL/SV/ALL, ax1: TP/FP/FN, ax2: QUAL - min_qual
    std::vector< std::vector< std::vector<float> > > query;

    /// truth counts; ax0: SNP/INDEL/SV/ALL, ax1: TP/FP/FN, ax2: QUAL - min_qual
    std::vector< std::vector< std::vector<float> > > truth;
};

/** @brief Tallies query and truth variant counts at each quality threshold, across all contigs. */
pr_counts tally_counts_by_qual(const std::unique_ptr<phaseblockData> & phasedata_ptr,
        int min_qual, int max_qual);

/**
 * @struct prec_recall_f1
 * @brief Precision, recall, and F1 score for one set of variant counts.
 */
struct prec_recall_f1 {
    float precision; ///< query true positive fraction; 1 if there are no query variants
    float recall;    ///< truth true positive fraction; 1 if there are no truth variants
    float f1;        ///< harmonic mean of precision and recall; 0 if their sum is not positive
};

/** @brief Computes precision, recall, and F1 score from query and truth variant counts. */
prec_recall_f1 compute_pr_f1(int query_tp, int query_fp, int truth_tp, int truth_fn);

/** @brief Writes precision-recall TSV results and prints console summary. */
void write_results(std::unique_ptr<phaseblockData> & phasedata_ptr);

/** @brief Writes all pipeline configuration parameters to TSV file. */
void write_params();

#endif
