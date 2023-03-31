#ifndef _VARIANT_H_
#define _VARIANT_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#include "htslib/vcf.h"

#include "fasta.h"
#include "defs.h"

class ctgVariants {
public:

    // constructor
    ctgVariants();

    // helper functions
    void add_cluster(int g);
    void add_var(int pos, int rlen, uint8_t hap, uint8_t type, uint8_t loc,
            std::string ref, std::string alt, uint8_t orig_gt,
            float gq, float vq);

    // data (all of size n)
    std::vector<int> poss;          // variant start positions (0-based)
    std::vector<int> rlens;         // reference lengths
    std::vector<uint8_t> haps;      // variant haplotype
    std::vector<uint8_t> types;     // variant type: NONE, SUB, INS, DEL, GRP
    std::vector<uint8_t> locs;      // BED location: INSIDE, OUTSIDE, BORDER
    std::vector<std::string> refs;  // variant reference allele
    std::vector<std::string> alts;  // variant alternate allele (always one)
    std::vector<uint8_t> orig_gts;  // simple genotype (0|1, 1|0, or 1|1)
    std::vector<float> gt_quals;    // genotype quality (0-60)
    std::vector<float> var_quals;   // variant quality (0-60)
    int n = 0;

    // set during (sw_)cluster()
    std::vector<int> clusters;      // indices of clusters in this struct's vectors

    // set during prec_recall_aln()
    std::vector<uint8_t> errtypes;  // error type: TP, FP, FN, PP
    std::vector<float> callq;       // call quality (for truth, of associated call)
                                    // or for call, min quality in sync group
    std::vector<float> credit;      // fraction of TP for partial positive (PP)
};

class variantData {
public:
    // constructors
    variantData();
    variantData(std::string vcf_fn, 
            std::shared_ptr<fastaData> reference, int callset);

    // functions
    void write_vcf(std::string vcf_fn);
    void print_variant(FILE* out_fp, std::string ctg, int pos, int type,
        std::string ref, std::string alt, float qual, std::string gt);
    void set_header(const std::unique_ptr<variantData> & vcf);
    void add_variants( const std::vector<int> & cigar, 
        int hap, int ref_pos,
        const std::string & ctg, 
        const std::string & query, 
        const std::string & ref,
        int qual);

    // data
    int callset;                     // 0=QUERY, 1=TRUTH
    std::string filename;
    std::string sample;
    std::vector<std::string> contigs;
    std::vector<int> lengths;
    std::vector< // ctg_variants[hap][ctg] -> variants
        std::unordered_map<
            std::string, 
            std::shared_ptr<ctgVariants> > > ctg_variants;

    std::shared_ptr<fastaData> ref;
};

#endif
