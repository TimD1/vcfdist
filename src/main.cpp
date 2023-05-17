#include "variant.h"
#include "print.h"
#include "globals.h"
#include "fasta.h"
#include "bed.h"
#include "dist.h"
#include "cluster.h"
#include "phase.h"
#include "timer.h"

Globals g;
std::vector<std::string> type_strs = {"REF", "SNP", "INS", "DEL", "CPX"};
std::vector<std::string> type_strs2 = {"ALL", "SNP", "INS", "DEL", "INDEL"};
std::vector<std::string> vartype_strs = {"SNP", "INDEL"};
std::vector<std::string> error_strs = {"TP", "FP", "FN", "PP", "PE", "GE", "??"};
std::vector<std::string> gt_strs = {
        "0", "1", "0|0", "0|1", "1|0", "1|1", "1|2", "2|1", ".|.", "M|N" };
std::vector<std::string> region_strs = {"OUTSIDE", "INSIDE ", "BORDER ", "OFF CTG"};
std::vector<std::string> aln_strs = {"QUERY1-TRUTH1", "QUERY1-TRUTH2", "QUERY2-TRUTH1", "QUERY2-TRUTH2"};
std::vector<std::string> callset_strs = {"QUERY", "TRUTH"};
std::vector<std::string> phase_strs = {"=", "X", "?"};
std::vector<std::string> timer_strs = {"reading", "clustering", "realigning",
    "reclustering", "superclustering", "aligning", "P/R", "S-W", "phasing", "writing", "total"};
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);
    g.init_timers(timer_strs);

    // parse reference fasta
g.timers[TIME_TOTAL].start();
g.timers[TIME_READ].start();
    std::shared_ptr<fastaData> ref_ptr(new fastaData(g.ref_fasta_fp));

    // parse query and truth VCFs
    std::unique_ptr<variantData> query_ptr(
            new variantData(g.query_vcf_fn, ref_ptr, QUERY));
    std::unique_ptr<variantData> truth_ptr(
            new variantData(g.truth_vcf_fn, ref_ptr, TRUTH));
g.timers[TIME_READ].stop();
g.timers[TIME_WRITE].start();
    query_ptr->write_vcf(g.out_prefix + "orig-query.vcf");
    truth_ptr->write_vcf(g.out_prefix + "orig-truth.vcf");
g.timers[TIME_WRITE].stop();

    // ensure each input contains all contigs in BED
    check_contigs(query_ptr, truth_ptr, ref_ptr);

    // cluster, realign, and cluster query VCF
    if (!g.keep_query) {
g.timers[TIME_CLUST].start();
        g.simple_cluster ? gap_cluster(query_ptr, QUERY) : swg_cluster(query_ptr, 
                g.query_sub, g.query_open, g.query_extend, QUERY); 
g.timers[TIME_CLUST].stop();
g.timers[TIME_REALN].start();
        query_ptr = wf_swg_realign(query_ptr, ref_ptr, 
                g.query_sub, g.query_open, g.query_extend, QUERY);
        query_ptr->left_shift();
g.timers[TIME_REALN].stop();
    }
g.timers[TIME_RECLUST].start();
    g.simple_cluster ? gap_cluster(query_ptr, QUERY) : swg_cluster(query_ptr, 
            g.query_sub, g.query_open, g.query_extend, QUERY); 
g.timers[TIME_RECLUST].stop();

    // cluster, realign, and cluster truth VCF
    if (!g.keep_truth) {
g.timers[TIME_CLUST].start();
        g.simple_cluster ?  gap_cluster(truth_ptr, TRUTH) : swg_cluster(truth_ptr, 
                g.truth_sub, g.truth_open, g.truth_extend, TRUTH); 
g.timers[TIME_CLUST].stop();
g.timers[TIME_REALN].start();
        truth_ptr = wf_swg_realign(truth_ptr, ref_ptr, 
                g.truth_sub, g.truth_open, g.truth_extend, TRUTH);
        truth_ptr->left_shift();
g.timers[TIME_REALN].stop();
    }
g.timers[TIME_RECLUST].start();
    g.simple_cluster ? gap_cluster(truth_ptr, TRUTH) : swg_cluster(truth_ptr, 
            g.truth_sub, g.truth_open, g.truth_extend, TRUTH); 
g.timers[TIME_RECLUST].stop();

    if (g.exit) {
g.timers[TIME_WRITE].start();
        query_ptr->write_vcf(g.out_prefix + "query.vcf");
        truth_ptr->write_vcf(g.out_prefix + "truth.vcf");
g.timers[TIME_WRITE].stop();
g.timers[TIME_TOTAL].stop();
        INFO(" ")
        INFO("TIMERS")
        for (int i = 0; i < TIME_TOTAL+1; i++) { g.timers[i].print(); }
        return EXIT_SUCCESS;
    }

    // calculate superclusters
g.timers[TIME_SUPCLUST].start();
    std::shared_ptr<superclusterData> clusterdata_ptr(
            new superclusterData(query_ptr, truth_ptr, ref_ptr));
g.timers[TIME_SUPCLUST].stop();

    // calculate precision/recall, edit distance, and local phasing
g.timers[TIME_ALIGN].start();
    editData edits = alignment_wrapper(clusterdata_ptr);
g.timers[TIME_ALIGN].stop();

    // calculate global phasings
g.timers[TIME_PHASE].start();
    std::unique_ptr<phaseData> phasedata_ptr(new phaseData(clusterdata_ptr));
g.timers[TIME_PHASE].stop();

    // write supercluster/phaseblock results in CSV format
g.timers[TIME_WRITE].start();
    write_results(phasedata_ptr, edits);

    // save new VCF
    query_ptr->write_vcf(g.out_prefix + "query.vcf");
    truth_ptr->write_vcf(g.out_prefix + "truth.vcf");
    phasedata_ptr->write_summary_vcf(g.out_prefix + "summary.vcf");
g.timers[TIME_WRITE].stop();

    // report timing results
g.timers[TIME_TOTAL].stop();
    INFO(" ")
    INFO("TIMERS")
    for (int i = 0; i < TIME_TOTAL+1; i++) { g.timers[i].print(); }
    return EXIT_SUCCESS;
}
