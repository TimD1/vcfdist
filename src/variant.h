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
    ctgVariants(const std::string & ctg);

    // helper functions
    void remove_vars(const std::vector<int> & indices);
    void add_var(std::shared_ptr<ctgVariants> other_vars, int idx);
    void add_var(int pos, int rlen, uint8_t type, uint8_t loc,
        std::string ref, std::string alt, uint8_t orig_gt, float gt_qual, float var_qual, 
        int phase_set, int supercluster = -1, uint8_t calc_gt = GT_REF_REF, 
        uint8_t hap1_errtype = ERRTYPE_UN, uint8_t hap2_errtype = ERRTYPE_UN, 
        int hap1_sync_group = 0, int hap2_sync_group = 0, float hap1_callq = 0, float hap2_callq = 0, 
        int hap1_ref_ed = 0, int hap2_ref_ed = 0, int hap1_query_ed = 0, int hap2_query_ed = 0, 
        float hap1_credit = 0, float hap2_credit = 0);
    void print_var_info(FILE* out_vcf, std::shared_ptr<fastaData> ref, 
            const std::string & ctg, int idx);
    void print_var_empty(FILE* out_vcf, int sc_idx, int phase_block, bool query = false);
    void print_var_sample(FILE* out_vcf, int var_idx, int hap_idx, const std::string & gt, int sc_idx, 
            int phase_block, bool phase_switch, bool phase_flip, bool query = false);
    bool var_on_hap(int var_idx, int hap_idx, bool calc = false) const;
    void set_var_calcgt_on_hap(int var_idx, int hap, bool set = true, bool ignore_errors = false);
    int get_vartype(int var_idx);
    int set_allele_errtype(int var_idx);
    bool calcgt_is_swapped(int var_idx) const;

    // originally parsed data (size n)
    std::string ctg;
    std::vector<int> poss;          // variant start positions (0-based)
    std::vector<int> rlens;         // reference lengths
    std::vector<uint8_t> types;     // variant type: NONE, SUB, INS, DEL, CPX
    std::vector<uint8_t> locs;      // BED location: INSIDE, OUTSIDE, BORDER
    std::vector<std::string> refs;  // variant reference allele
    std::vector<std::string> alts;  // variant alternate allele (always one)
    std::vector<uint8_t> orig_gts;  // simple genotype (0|1, 1|0, or 1|1)
    std::vector<float> gt_quals;    // genotype quality (0-60)
    std::vector<float> var_quals;   // variant quality (0-60)
    std::vector<int> phase_sets;    // integer representing variant phase set (0 = missing)
    std::vector<int> superclusters; // initially -1, set during superclustering
    int n = 0;

    // set during clustering (size nc+1)
    // TODO: these are now ONLY used during clustering, should be removed from the class
    std::vector<int> clusters;      // variant indices of clusters in this struct's vectors
    std::vector<int> left_reaches;  // cluster leftmost reach
    std::vector<int> right_reaches; // cluster rightmost reach
    int nc = 0;

    // set during prec_recall_aln() (size (2, n), additional axis for haplotype)
    std::vector<uint8_t> calc_gts;  // calculated genotype (0|1, 1|0, or 1|1), only set for query 
    std::vector< std::vector<uint8_t> > errtypes;  // error type: TP, FP, FN
    std::vector< std::vector<int> > sync_group;    // group of variants that participate in credit
    std::vector< std::vector<float> > callq;   // min call quality in sync group (for truth, of associated call)
    std::vector< std::vector<int> > ref_ed;    // reference edit distance in sync group
    std::vector< std::vector<int> > query_ed;  // query edit distance in sync group
    std::vector< std::vector<float> > credit;  // percentage reduction in edit dist (ref->query)

    // set during phase() (size n)
    std::vector<int> phases;             // variant keep/swap/unknown, from alignment (calc_gt relative to orig_gt)
    std::vector<int> pb_phases;          // phaseblock keep/swap, from phasing algorithm
    std::vector<int> ac_errtype;         // allele count error type (e.g. 0|1 -> 1|1)
};

class variantData {
public:
    // constructors
    variantData();

    // functions
    void write_vcf(std::string vcf_fn);
    void print_variant(FILE* out_fp, const std::string & ctg, int pos, int type,
        const std::string & ref, const std::string & alt, float qual, int phase_set,
        const std::string & gt);
    void set_header(const std::shared_ptr<variantData> vcf);
    void add_variants(const std::vector<int> & cigar, int hap,
            int ref_pos, const std::string & ctg, const std::string & query, 
            const std::string & ref, int qual, int phase_set);
    void merge(std::shared_ptr<variantData> other_variants);
    void print_phase_info(int callset);

    // data
    std::shared_ptr<fastaData> ref;
    int callset;                     // 0=QUERY, 1=TRUTH
    std::string filename;

    std::string sample;
    std::vector<std::string> contigs;
    std::vector<int> lengths;
    std::vector<int> ploidy;
    std::vector< // variants[hap][ctg] -> ctg_variants
        std::unordered_map<
            std::string, 
            std::shared_ptr<ctgVariants> > > variants;
};

void parse_variants(const std::string & vcf_fn, 
        std::shared_ptr<variantData> variant_data,
        std::shared_ptr<variantData> large_variant_data,
        std::shared_ptr<fastaData> reference, int callset);

#endif
