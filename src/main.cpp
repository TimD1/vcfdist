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
    std::shared_ptr<variantData> query_ptr(
            new variantData(g.query_vcf_fn, ref_ptr, QUERY));
    query_ptr->print_phase_info(QUERY);
    std::shared_ptr<variantData> truth_ptr(
            new variantData(g.truth_vcf_fn, ref_ptr, TRUTH));
    g.timers[TIME_READ].stop();
    truth_ptr->print_phase_info(TRUTH);

    // ensure each input contains all contigs in BED
    check_contigs(query_ptr, truth_ptr, ref_ptr);

    // cluster query VCF
    g.timers[TIME_CLUST].start();
    if (g.cluster_method == "gap" || g.cluster_method == "size") {
        simple_cluster(query_ptr, QUERY);
    } else if (g.cluster_method == "biwfa") {
        if (g.verbosity >= 1) INFO(" ");
        if (g.verbosity >= 1) INFO("%s[Q %d/%d] Wavefront clustering %s VCF%s '%s'", 
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
    } else { 
        ERROR("Unexpected clustering method '%s'", g.cluster_method.data()); 
    }
    g.timers[TIME_CLUST].stop();

    // cluster truth VCF
    g.timers[TIME_CLUST].start();
    if (g.cluster_method == "gap" || g.cluster_method == "size") {
        simple_cluster(truth_ptr, TRUTH); 
    } else if (g.cluster_method == "biwfa") {
        if (g.verbosity >= 1) INFO(" ");
        if (g.verbosity >= 1) INFO("%s[T %d/%d] Wavefront clustering %s VCF%s '%s'", 
                COLOR_PURPLE, TIME_CLUST, TIME_TOTAL-1, callset_strs[TRUTH].data(), 
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
