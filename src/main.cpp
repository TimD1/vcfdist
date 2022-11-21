#include "variant.h"
#include "print.h"
#include "globals.h"
#include "fasta.h"
#include "bed.h"
#include "dist.h"
#include "cluster.h"
#include "phase.h"

Globals g;
std::vector<std::string> type_strs = {"REF", "SUB", "INS", "DEL", "GRP"};
std::vector<std::string> vartype_strs = {"SNP", "INDEL"};
std::vector<std::string> error_strs = {"TP", "FP", "FN", "PP", "PE", "GE", "??"};
std::vector<std::string> gt_strs = {
".|.", "0|0", "0|1", "0|2", "1|0", "1|1", "1|2", "2|0", "2|1", "2|2", "?|?" };
std::vector<std::string> region_strs = {"OUTSIDE", "INSIDE ", "BORDER ", "OFF CTG"};
std::vector<std::string> aln_strs = {"QUERY1-TRUTH1", "QUERY1-TRUTH2", "QUERY2-TRUTH1", "QUERY2-TRUTH2"};
std::vector<std::string> phase_strs = {"=", "X", "?"};
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);

    // parse reference fasta
    std::shared_ptr<fastaData> ref_ptr(new fastaData(g.ref_fasta_fp));

    // parse, realign, and cluster query VCF
    std::unique_ptr<variantData> query_ptr(
            new variantData(g.query_vcf_fn, ref_ptr));
    query_ptr->write_vcf(g.out_prefix + "orig_query.vcf");
    g.simple_cluster ? cluster(query_ptr) :
        sw_cluster(query_ptr, g.query_sub, g.query_open, g.query_extend); 
    if (!g.keep_query) {
        query_ptr = sw_realign(query_ptr, ref_ptr, 
                g.query_sub, g.query_open, g.query_extend);
        g.simple_cluster ?  cluster(query_ptr) :
            sw_cluster(query_ptr, g.query_sub, g.query_open, g.query_extend); 
    }
    
    // parse, realign, and cluster truth VCF
    std::unique_ptr<variantData> truth_ptr(
            new variantData(g.truth_vcf_fn, ref_ptr));
    truth_ptr->write_vcf(g.out_prefix + "orig_truth.vcf");
    g.simple_cluster ?  cluster(truth_ptr) :
        sw_cluster(truth_ptr, g.truth_sub, g.truth_open, g.truth_extend); 
    if (!g.keep_truth) {
        truth_ptr = sw_realign(truth_ptr, ref_ptr, 
                g.truth_sub, g.truth_open, g.truth_extend);
        g.simple_cluster ?  cluster(truth_ptr) :
            sw_cluster(truth_ptr, g.truth_sub, g.truth_open, g.truth_extend); 
    }

    if (g.exit) {
        query_ptr->write_vcf(g.out_prefix + "query.vcf");
        truth_ptr->write_vcf(g.out_prefix + "truth.vcf");
        return EXIT_SUCCESS;
    }

    // calculate superclusters
    std::shared_ptr<superclusterData> clusterdata_ptr(
            new superclusterData(query_ptr, truth_ptr, ref_ptr));

    // calculate edit distance and local phasing
    int dist = alignment_wrapper(clusterdata_ptr);

    // calculate global phasings
    std::unique_ptr<phaseData> phasedata_ptr(new phaseData(clusterdata_ptr));

    // write supercluster/phaseblock results in CSV format
    write_results(phasedata_ptr, dist);

    // save new VCF
    query_ptr->write_vcf(g.out_prefix + "query.vcf");
    truth_ptr->write_vcf(g.out_prefix + "truth.vcf");

    return EXIT_SUCCESS;
}
