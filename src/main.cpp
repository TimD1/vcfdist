#include "vcf.h"
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
std::vector<std::string> aln_strs = {"CAL1-HAP1", "CAL1-HAP2", "CAL2-HAP1", "CAL2-HAP2"};
std::vector<std::string> phase_strs = {"0", "1", "-"};
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);

    // parse reference fasta
    fastaData ref(g.ref_fasta_fp);

    // parse calls VCF, cluster variants, get min edit dist
    vcfData calls(g.calls_vcf_fn, &ref);
    cluster(&calls);
    vcfData calls_min_ed = edit_dist_realign(&calls, &ref);
    calls.write(g.out_prefix + "orig_calls.vcf");
    calls_min_ed.write(g.out_prefix + "calls.vcf");

    // parse ground truth VCF, cluster variants, get min edit dist
    vcfData truth(g.truth_vcf_fn, &ref);
    cluster(&truth);
    vcfData truth_min_ed = edit_dist_realign(&truth, &ref);
    truth.write(g.out_prefix + "orig_truth.vcf");
    truth_min_ed.write(g.out_prefix + "truth.vcf");

    // save per-cluster alignment info
    clusterData clusters = edit_dist(&calls, &truth, &ref);

    // phase clusters
    phaseData phasings(&clusters);

    // store results in CSV format
    write_results(phasings);

    return EXIT_SUCCESS;
}
