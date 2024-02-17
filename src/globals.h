#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "bed.h"
#include "defs.h"
#include "timer.h"

class Globals {
public:
    // constructors
    Globals() {;}

    // input files
    std::string ref_fasta_fn;
    FILE* ref_fasta_fp;
    std::string query_vcf_fn;
    std::string truth_vcf_fn;
    std::string bed_fn;
    bedData bed;
    bool bed_exists = false;
    bool write = true;

    // variant filtering
    std::vector<std::string> filters;
    std::vector<int> filter_ids;
    int min_qual = 0;
    int max_qual = 60;
    int min_size = 1;
    int max_size = 5000;

    // clustering
    std::string cluster_method = "biwfa";
    int cluster_min_gap = 50;
    int reach_min_gap = 10;
    int max_cluster_itrs = 4;

    // re-alignment
    bool realign_truth = false;
    bool realign_query = false;
    bool realign_only = false;
    int sub = 5;
    int open = 6;
    int extend = 2;

    // phasing
    double phase_threshold = 0.6;

    // precision-recall
    double credit_threshold = 0.7;

    // edit distance
    bool distance = false;
    int eval_sub = 3;
    int eval_open = 2;
    int eval_extend = 1;

    // memory params
    int max_threads = 64;
    double max_ram = 64; // GB
    int thread_nsteps;
    std::vector<int> thread_steps;
    std::vector<float> ram_steps;

    // high-level options
    bool advanced = false;
    int verbosity = 1;
    int sv_threshold = 50;
    std::string out_prefix;
    std::string cmd;

    // member functions
    void parse_args(int argc, char ** argv);
    void print_version() const;
    void print_usage() const; 
    void print_citation() const; 
    void init_timers(std::vector<std::string> timer_strs);

    // program data
    const std::string VERSION = "2.3.4";
    const std::string PROGRAM = "vcfdist";
    std::vector<timer> timers;
};

std::string parent_path(std::string out_prefix);
void create_directory(std::string dir);

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
extern std::vector<std::string> switch_strs;

#endif
