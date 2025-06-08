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
    std::shared_ptr<variantData> query_ptr(new variantData());
    parse_variants(g.query_vcf_fn, query_ptr, ref_ptr, QUERY);

    std::shared_ptr<variantData> truth_ptr(new variantData());
    parse_variants(g.truth_vcf_fn, truth_ptr, ref_ptr, TRUTH);
    g.timers[TIME_READ].stop();

    // write parsed per-variant information
    g.timers[TIME_WRITE].start();
    if (g.write) {
        query_ptr->write_vcf(g.out_prefix + "query.vcf");
        truth_ptr->write_vcf(g.out_prefix + "truth.vcf");
    }
    g.timers[TIME_WRITE].stop();

    // ensure each input contains all contigs in BED
    check_contigs(query_ptr, truth_ptr, ref_ptr);

    // cluster query VCF
    g.timers[TIME_CLUSTER].start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[Q %d/%d] Clustering %s VCF%s '%s'", 
            COLOR_PURPLE, TIME_CLUSTER, TIME_TOTAL-1, callset_strs[QUERY].data(), 
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

    // cluster truth VCF
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[T %d/%d] Clustering %s VCF%s '%s'", 
            COLOR_PURPLE, TIME_CLUSTER, TIME_TOTAL-1, callset_strs[TRUTH].data(), 
            COLOR_WHITE, truth_ptr->filename.data());
    threads.clear();
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

    // superclustering: merge per-hap variant info in constructor, then supercluster
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%d/%d] Superclustering TRUTH and QUERY variants%s",
            COLOR_PURPLE, TIME_CLUSTER, TIME_TOTAL-1, COLOR_WHITE);
    std::shared_ptr<superclusterData> sc_data_ptr(
            new superclusterData(query_ptr, truth_ptr, ref_ptr));
    sc_data_ptr->supercluster(false);
    auto sc_groups = sort_superclusters(sc_data_ptr);
    g.timers[TIME_CLUSTER].stop();

    // evaluation: precision/recall and genotypes
    g.timers[TIME_ALIGN_EVAL].start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%d/%d] Evaluating variant calls %s",
            COLOR_PURPLE, TIME_ALIGN_EVAL, TIME_TOTAL-1, COLOR_WHITE);
    precision_recall_threads_wrapper(sc_data_ptr, sc_groups);
    INFO("    done with variant call evaluation");
    g.timers[TIME_ALIGN_EVAL].stop();

    // calculate phasing statistics
    g.timers[TIME_PHASE].start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%d/%d] Phasing superclusters%s",
            COLOR_PURPLE, TIME_PHASE, TIME_TOTAL-1, COLOR_WHITE);
    std::unique_ptr<phaseblockData> phasedata_ptr(new phaseblockData(sc_data_ptr));
    g.timers[TIME_PHASE].stop();

    // write phasing results
    g.timers[TIME_WRITE].start();
    if (g.write) phasedata_ptr->write_switchflips();
    write_results(phasedata_ptr);
    if (g.write) phasedata_ptr->write_summary_vcf(g.out_prefix + "summary.vcf");
    g.timers[TIME_WRITE].stop();

    // report timing results
    g.timers[TIME_TOTAL].stop();
    write_runtime();
    if (g.verbosity >= 1) {
        INFO(" ")
        INFO("Timers:")
        for (int i = 0; i <= TIME_TOTAL; i++) { g.timers[i].print(i); }
    }
    return EXIT_SUCCESS;
}
