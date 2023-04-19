#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <filesystem>

#include "bed.h"
#include "defs.h"

class Globals {
public:
    // command-line data, set using parse_args()
    std::string ref_fasta_fn;
    FILE* ref_fasta_fp;

    std::string query_vcf_fn;
    std::string truth_vcf_fn;

    std::string bed_fn;
    bedData bed;
    bool bed_exists = false;

    std::string out_prefix;
    std::string cmd;

    int cluster_min_gap = 50;
    int max_cluster_itrs = 4;
    int verbosity = 1;
    int min_qual = 0;
    int max_qual = 60;
    int min_size = 1;
    int max_size = 5000;

    bool exit = false;
    bool advanced = false;

    int query_sub = 5;
    int query_open = 6;
    int query_extend = 2;
    int truth_sub = 5;
    int truth_open = 6;
    int truth_extend = 2;
    int eval_sub = 3;
    int eval_open = 2;
    int eval_extend = 1;
    std::vector<bool> query_penalties_set = {false, false, false};
    std::vector<bool> truth_penalties_set = {false, false, false};
    std::vector<bool> eval_penalties_set  = {false, false, false};
    
    bool keep_truth = false;
    bool keep_query = false;
    bool simple_cluster = false;

    // constructors
    Globals() {;}

    // member functions
    void parse_args(int argc, char ** argv);
    void print_version() const;
    void print_usage() const; 

    // program data
    const std::string VERSION = "1.2.2";
    const std::string PROGRAM = "vcfdist";
};

extern Globals g;

// defined in main.cpp
extern std::vector<std::string> callset_strs;
extern std::vector<std::string> vartype_strs;
extern std::vector<std::string> type_strs;
extern std::vector<std::string> type_strs2;
extern std::vector<std::string> error_strs;
extern std::vector<std::string> gt_strs;
extern std::vector<std::string> region_strs;
extern std::vector<std::string> aln_strs;
extern std::vector<std::string> phase_strs;

#endif
