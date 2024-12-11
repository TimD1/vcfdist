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
    std::shared_ptr<variantData> large_query_ptr(new variantData());
    parse_variants(g.query_vcf_fn, query_ptr, large_query_ptr, ref_ptr, QUERY);
    query_ptr->print_phase_info(QUERY);
    std::shared_ptr<variantData> truth_ptr(new variantData());
    std::shared_ptr<variantData> large_truth_ptr(new variantData());
    parse_variants(g.truth_vcf_fn, truth_ptr, large_truth_ptr, ref_ptr, TRUTH);
    truth_ptr->print_phase_info(TRUTH);
    g.timers[TIME_READ].stop();

    // ensure each input contains all contigs in BED
    check_contigs(query_ptr, truth_ptr, ref_ptr);

    // cluster query VCF
    g.timers[TIME_CLUST].start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[Q %d/%d] Clustering %s VCF%s '%s'", 
            COLOR_PURPLE, TIME_CLUST, TIME_TOTAL-1, callset_strs[QUERY].data(), 
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
    g.timers[TIME_CLUST].stop();

    // cluster truth VCF
    g.timers[TIME_CLUST].start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[T %d/%d] Clustering %s VCF%s '%s'", 
            COLOR_PURPLE, TIME_CLUST, TIME_TOTAL-1, callset_strs[TRUTH].data(), 
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
    g.timers[TIME_CLUST].stop();

    // calculate superclusters
    g.timers[TIME_SUPCLUST].start();
    std::shared_ptr<superclusterData> clusterdata_ptr(
            new superclusterData(query_ptr, truth_ptr, ref_ptr));

    // calculate supercluster sizes
    auto sc_groups = sort_superclusters(clusterdata_ptr);
    g.timers[TIME_SUPCLUST].stop();

    // calculate precision/recall and genotypes
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
        query_ptr->write_vcf(g.out_prefix + "query.vcf");
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
