#include "variant.h"
#include "print.h"
#include "globals.h"
#include "fasta.h"
#include "bed.h"
#include "dist.h"
#include "cluster.h"
#include "phase.h"

Globals g;
std::vector<std::string> type_strs = { "REF", "SUB", "INS", "DEL", "GRP"};
std::vector<std::string> error_strs = { "TP", "FP", "FN", "PP", "PE", "GE"};
std::vector<std::string> gt_strs = {
".|.", "0|0", "0|1", "0|2", "1|0", "1|1", "1|2", "2|0", "2|1", "2|2", "?|?" };
std::vector<std::string> region_strs = {"OUTSIDE", "INSIDE ", "BORDER ", "OFF CTG"};
std::vector<std::string> aln_strs = {"CALLS1-TRUTH1", "CALLS1-TRUTH2", "CALLS2-TRUTH1", "CALLS2-TRUTH2"};
std::vector<std::string> phase_strs = {"=", "X", "?"};
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);

    // parse reference fasta
    std::shared_ptr<fastaData> ref_ptr(new fastaData(g.ref_fasta_fp));

    // parse calls VCF, cluster variants, get min edit dist
    std::unique_ptr<variantData> calls_ptr(new variantData(g.calls_vcf_fn, ref_ptr));
    cluster(calls_ptr);
    variantData calls_min_ed = edit_dist_realign(calls_ptr, ref_ptr);
    calls_ptr->write(g.out_prefix + "orig_calls.vcf");
    calls_min_ed.write(g.out_prefix + "calls.vcf");

    // parse ground truth VCF, cluster variants, get min edit dist
    std::unique_ptr<variantData> truth_ptr(new variantData(g.truth_vcf_fn, ref_ptr));
    cluster(truth_ptr);
    variantData truth_min_ed = edit_dist_realign(truth_ptr, ref_ptr);
    truth_ptr->write(g.out_prefix + "orig_truth.vcf");
    truth_min_ed.write(g.out_prefix + "truth.vcf");

    // calculate superclusters
    std::shared_ptr<clusterData> clusterdata_ptr(new clusterData(calls_ptr, truth_ptr, ref_ptr));

    // calculate edit distance and local phasing
    alignment_wrapper(clusterdata_ptr);

    // calculate global phasings
    std::unique_ptr<phaseData> phasedata_ptr(new phaseData(clusterdata_ptr));

    // store results in CSV format
    write_results(phasedata_ptr);

    return EXIT_SUCCESS;
}
