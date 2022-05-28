#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "htslib/vcf.h"

#include "bed.h"

class Globals {
public:
    // command-line data, set using parse_args()
    std::string ref_fasta_fn;
    FILE* ref_fasta_fp;

    std::string calls_vcf_fn;
    htsFile * calls_vcf_fp;

    std::string truth_vcf_fn;
    htsFile * truth_vcf_fp;

    std::string bed_fn;
    bedData bed;

    std::string out_dir;

    int gap = 50;
    int print_verbosity = 0;

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
extern std::vector< std::string > type_strs;
extern std::vector< std::string > gt_strs;
extern std::vector< std::string > region_strs;

#endif
