#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <filesystem>

#include "bed.h"

#define PEN_SUB    0 // penalties order in boolean array *_penalties_set
#define PEN_OPEN   1
#define PEN_EXTEND 2

class Globals {
public:
    // command-line data, set using parse_args()
    std::string ref_fasta_fn;
    FILE* ref_fasta_fp;

    std::string calls_vcf_fn;
    std::filesystem::path calls_vcf_path;

    std::string truth_vcf_fn;
    std::filesystem::path truth_vcf_path;

    std::string bed_fn;
    bedData bed;

    std::string out_prefix;

    int cluster_min_gap = 50;
    int print_verbosity = 0;
    int min_qual = 0;

    int calls_sub = 1;
    int calls_open = 1;
    int calls_extend = 1;
    int truth_sub = 1;
    int truth_open = 1;
    int truth_extend = 1;
    std::vector<bool> calls_penalties_set = {false, false, false};
    std::vector<bool> truth_penalties_set = {false, false, false};
    
    bool keep_truth = false;
    bool keep_calls = false;
    bool simple_cluster = false;

    // constructors
    Globals() {;}

    // member functions
    void parse_args(int argc, char ** argv);
    void print_version() const;
    void print_usage() const; 

    // program data
    const std::string VERSION = "0.0.1";
    const std::string PROGRAM = "vcfdist";
};

extern Globals g;
extern std::vector<std::string> type_strs;
extern std::vector<std::string> error_strs;
extern std::vector<std::string> gt_strs;
extern std::vector<std::string> region_strs;
extern std::vector<std::string> aln_strs;
extern std::vector<std::string> phase_strs;

#endif
