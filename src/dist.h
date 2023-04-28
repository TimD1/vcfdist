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
        std::unordered_map<idx1, int> & aln_ptrs, 
        std::vector< std::unordered_map<idx1, idx1> > & pred_map,
        std::vector<int> & pr_query_ref_end, bool print
        );

void get_prec_recall_path_sync(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::unordered_map<idx1, int> & aln_ptrs, 
        std::unordered_map<idx1, int> & path_ptrs, 
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
        const std::unordered_map<idx1, int> & aln_ptrs, 
        const std::unordered_map<idx1, int> & path_ptrs, 
        const std::unordered_map<idx1, int> & path_scores, 
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
        const std::string & ctg, 
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        const std::vector< std::vector<idx1> > & path,
        const std::vector< std::vector<bool> > & sync,
        const std::vector< std::vector<bool> > & edits,
        const std::unordered_map<idx1, int> & aln_ptrs, 
        const std::unordered_map<idx1, int> & path_ptrs, 
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
