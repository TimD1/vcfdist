#include <thread>

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
std::vector<std::string> vartype_strs = {"SNP", "INDEL", "SV", "ALL"};
std::vector<std::string> error_strs = {"TP", "FP", "FN", "PE", "GE", "??"};
std::vector<std::string> gt_strs = {
        "0", "1", "0|0", "0|1", "1|0", "1|1", "1|2", "2|1", ".|.", "M|N" };
std::vector<std::string> region_strs = {"OUTSIDE", "INSIDE ", "BORDER ", "OFF CTG"};
std::vector<std::string> aln_strs = {"QUERY1-TRUTH1", "QUERY1-TRUTH2", "QUERY2-TRUTH1", "QUERY2-TRUTH2"};
std::vector<std::string> callset_strs = {"QUERY", "TRUTH"};
std::vector<std::string> phase_strs = {"=", "X", "?"};
std::vector<std::string> timer_strs = {"reading", "clustering", "realigning", 
    "reclustering", "superclustering", "precision/recall", "phasing", "writing", "total"};
std::vector<std::string> switch_strs = {"FLIP", "SWITCH", "SWITCH+FLIP", "SWITCH_ERR", "FLIP_BEG", "FLIP_END", "NONE"};
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);
    g.init_timers(timer_strs);

    g.timers[TIME_WRITE].start();
    write_params();
    g.timers[TIME_WRITE].stop();

    // parse reference fasta
    g.timers[TIME_TOTAL].start();
    g.timers[TIME_READ].start();
    std::shared_ptr<fastaData> ref_ptr(new fastaData(g.ref_fasta_fp));

    // parse query and truth VCFs
    std::shared_ptr<variantData> query_ptr(
            new variantData(g.query_vcf_fn, ref_ptr, QUERY));
    query_ptr->print_phase_info(QUERY);
    std::shared_ptr<variantData> truth_ptr(
            new variantData(g.truth_vcf_fn, ref_ptr, TRUTH));
    g.timers[TIME_READ].stop();
    truth_ptr->print_phase_info(TRUTH);

    // write results
    if (g.write) {
        if (g.realign_query) {
            g.timers[TIME_WRITE].start();
            if (g.verbosity >= 1) INFO("  Writing original query VCF to '%s'", 
                    std::string(g.out_prefix + "orig-query.vcf").data());
            query_ptr->write_vcf(g.out_prefix + "orig-query.vcf");
            g.timers[TIME_WRITE].stop();
        }
        if (g.realign_truth) {
            g.timers[TIME_WRITE].start();
            if (g.verbosity >= 1) INFO("  Writing original truth VCF to '%s'", 
                    std::string(g.out_prefix + "orig-truth.vcf").data());
            truth_ptr->write_vcf(g.out_prefix + "orig-truth.vcf");
            g.timers[TIME_WRITE].stop();
        }
    }

    // ensure each input contains all contigs in BED
    check_contigs(query_ptr, truth_ptr, ref_ptr);

    // cluster, realign, and cluster query VCF
    if (g.realign_query) {
        g.timers[TIME_CLUST].start();
        if (g.cluster_method == "gap" || g.cluster_method == "size") {
            simple_cluster(query_ptr, QUERY);
        } else if (g.cluster_method == "biwfa") {
            if (g.verbosity >= 1) INFO(" ");
            if (g.verbosity >= 1) INFO("%s[Q 1/7] Wavefront clustering %s VCF%s '%s'", 
                    COLOR_PURPLE, callset_strs[QUERY].data(), 
                    COLOR_WHITE, query_ptr->filename.data());
            std::vector<std::thread> threads;
            for (int t = 0; t < HAPS*int(query_ptr->contigs.size()); t++) {
                threads.push_back(std::thread( wf_swg_cluster, 
                            query_ptr.get(), t/2 /* contig */, t%2, /* hap */
                            g.sub, g.open, g.extend)); 
                if ((t+1) % g.max_threads == 0) { // wait for thread batch to complete
                    for (std::thread & thread : threads) thread.join();
                    threads.clear();
                }
            }
            for (std::thread & thread : threads) thread.join();
        } else {
            ERROR("Unexpected clustering method '%s'", g.cluster_method.data());
        }
        g.timers[TIME_CLUST].stop();

        g.timers[TIME_REALN].start();
        query_ptr = wf_swg_realign(query_ptr, ref_ptr, 
                g.sub, g.open, g.extend, QUERY);
        query_ptr->left_shift();
        g.timers[TIME_REALN].stop();
    }

    if (g.realign_only && g.realign_query && g.write) { // realign only, exit early
        g.timers[TIME_WRITE].start();
        query_ptr->write_vcf(g.out_prefix + "query.vcf");
        g.timers[TIME_WRITE].stop();

    } else { // re-cluster based on new alignments
        g.timers[TIME_RECLUST].start();
        if (g.cluster_method == "gap" || g.cluster_method == "size") {
            simple_cluster(query_ptr, QUERY);
        } else if (g.cluster_method == "biwfa") {
            if (g.verbosity >= 1) INFO(" ");
            if (g.verbosity >= 1) INFO("%s[Q 3/7] Wavefront %sclustering %s VCF%s '%s'", 
                    COLOR_PURPLE, g.realign_query ? "re":"", callset_strs[QUERY].data(), 
                    COLOR_WHITE, query_ptr->filename.data());
            std::vector<std::thread> threads;
            for (int t = 0; t < HAPS*int(query_ptr->contigs.size()); t++) {
                threads.push_back(std::thread( wf_swg_cluster, 
                            query_ptr.get(), t/2 /* contig */, t%2, /* hap */
                            g.sub, g.open, g.extend)); 
                if ((t+1) % g.max_threads == 0) { // wait for thread batch to complete
                    for (std::thread & thread : threads) thread.join();
                    threads.clear();
                }
            }
            for (std::thread & thread : threads) thread.join();
        } else { 
            ERROR("Unexpected clustering method '%s'", g.cluster_method.data()); 
        }
        g.timers[TIME_RECLUST].stop();
    }

    // cluster, realign, and cluster truth VCF
    if (g.realign_truth) {
        g.timers[TIME_CLUST].start();
        if (g.cluster_method == "gap" || g.cluster_method == "size") {
            simple_cluster(truth_ptr, TRUTH);
        } else if (g.cluster_method == "biwfa") {
            if (g.verbosity >= 1) INFO(" ");
            if (g.verbosity >= 1) INFO("%s[T 1/7] Wavefront clustering %s VCF%s '%s'", 
                    COLOR_PURPLE, callset_strs[TRUTH].data(), 
                    COLOR_WHITE, truth_ptr->filename.data());
            std::vector<std::thread> threads;
            for (int t = 0; t < HAPS*int(truth_ptr->contigs.size()); t++) {
                threads.push_back(std::thread( wf_swg_cluster, 
                            truth_ptr.get(), t/2 /* contig */, t%2, /* hap */
                            g.sub, g.open, g.extend)); 
                if ((t+1) % g.max_threads == 0) { // wait for thread batch to complete
                    for (std::thread & thread : threads) thread.join();
                    threads.clear();
                }
            }
            for (std::thread & thread : threads) thread.join();
        } else { 
            ERROR("Unexpected clustering method '%s'", g.cluster_method.data()); 
        }
        g.timers[TIME_CLUST].stop();

        g.timers[TIME_REALN].start();
        truth_ptr = wf_swg_realign(truth_ptr, ref_ptr, 
                g.sub, g.open, g.extend, TRUTH);
        truth_ptr->left_shift();
        g.timers[TIME_REALN].stop();
    }

    if (g.realign_only) { // realign only, exit early
        if (g.realign_truth && g.write) {
            g.timers[TIME_WRITE].start();
            truth_ptr->write_vcf(g.out_prefix + "truth.vcf");
            g.timers[TIME_WRITE].stop();
        }
        g.timers[TIME_TOTAL].stop();
        if (g.verbosity >= 1) {
            INFO(" ")
            INFO("Timers:")
            for (int i = 0; i < TIME_TOTAL+1; i++) { g.timers[i].print(i); }
        }
        return EXIT_SUCCESS;

    } else {

        g.timers[TIME_RECLUST].start();
        if (g.cluster_method == "gap" || g.cluster_method == "size") {
            simple_cluster(truth_ptr, TRUTH); 
        } else if (g.cluster_method == "biwfa") {
            if (g.verbosity >= 1) INFO(" ");
            if (g.verbosity >= 1) INFO("%s[T 3/7] Wavefront %sclustering %s VCF%s '%s'", 
                    COLOR_PURPLE, g.realign_truth ? "re":"", callset_strs[TRUTH].data(), 
                    COLOR_WHITE, truth_ptr->filename.data());
            std::vector<std::thread> threads;
            for (int t = 0; t < HAPS*int(truth_ptr->contigs.size()); t++) {
                threads.push_back(std::thread( wf_swg_cluster, 
                            truth_ptr.get(), t/2 /* contig */, t%2, /* hap */
                            g.sub, g.open, g.extend)); 
                if ((t+1) % g.max_threads == 0) { // wait for thread batch to complete
                    for (std::thread & thread : threads) thread.join();
                    threads.clear();
                }
            }
            for (std::thread & thread : threads) thread.join();
        } else { 
            ERROR("Unexpected clustering method '%s'", g.cluster_method.data()); 
        }
        g.timers[TIME_RECLUST].stop();
    }

    // calculate superclusters
    g.timers[TIME_SUPCLUST].start();
    std::shared_ptr<superclusterData> clusterdata_ptr(
            new superclusterData(query_ptr, truth_ptr, ref_ptr));

    // calculate supercluster sizes
    auto sc_groups = sort_superclusters(clusterdata_ptr);
    g.timers[TIME_SUPCLUST].stop();

    // calculate precision/recall and local phasing
    g.timers[TIME_PR_ALN].start();
    precision_recall_threads_wrapper(clusterdata_ptr, sc_groups);
    INFO("    done with precision-recall");
    g.timers[TIME_PR_ALN].stop();

    // calculate global phasings
    g.timers[TIME_PHASE].start();
    std::unique_ptr<phaseblockData> phasedata_ptr(new phaseblockData(clusterdata_ptr));
    g.timers[TIME_PHASE].stop();

    // write supercluster/phaseblock results in CSV format
    g.timers[TIME_WRITE].start();
    if (g.write) phasedata_ptr->write_switchflips();
    write_results(phasedata_ptr);

    // save new VCF
    if (g.write) {
        if (g.realign_query)
            query_ptr->write_vcf(g.out_prefix + "query.vcf");
        if (g.realign_truth)
            truth_ptr->write_vcf(g.out_prefix + "truth.vcf");
        phasedata_ptr->write_summary_vcf(g.out_prefix + "summary.vcf");
    }
    g.timers[TIME_WRITE].stop();

    // report timing results
    g.timers[TIME_TOTAL].stop();
    if (g.verbosity >= 1) {
        INFO(" ")
        INFO("Timers:")
        for (int i = 0; i < TIME_TOTAL+1; i++) { g.timers[i].print(i); }
    }
    return EXIT_SUCCESS;
}
