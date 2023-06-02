#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <filesystem>

#include "bed.h"
#include "defs.h"
#include "timer.h"

class Globals {
public:
    // input files
    std::string ref_fasta_fn;
    FILE* ref_fasta_fp;
    std::string query_vcf_fn;
    std::string truth_vcf_fn;
    std::string bed_fn;
    bedData bed;
    bool bed_exists = false;

    // variant params
    int min_qual = 0;
    int max_qual = 60;
    int min_size = 1;
    int max_size = 5000;

    // clustering params
    bool keep_truth = false;
    bool keep_query = false;
    bool simple_cluster = false;
    int cluster_min_gap = 50;
    int reach_min_gap = 10;
    int max_cluster_itrs = 4;

    // memory params
    int max_threads = 64;
    float max_ram = 32; // GB
    int thread_nsteps;
    std::vector<int> thread_steps;
    std::vector<float> ram_steps;

    // high-level options
    bool exit = false;
    bool advanced = false;
    int verbosity = 1;
    std::string out_prefix;
    std::string cmd;

    // alignment parameters
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

    std::vector<timer> timers;

    // constructors
    Globals() {;}

    // member functions
    void parse_args(int argc, char ** argv);
    void print_version() const;
    void print_usage() const; 
    void print_citation() const; 
    void init_timers(std::vector<std::string> timer_strs);

    // program data
    const std::string VERSION = "1.3.0";
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
