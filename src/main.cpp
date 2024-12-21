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
    std::shared_ptr<variantData> query_ptr2(new variantData()); // second round eval
    parse_variants(g.query_vcf_fn, query_ptr, query_ptr2, ref_ptr, QUERY);
    query_ptr->print_phase_info(QUERY);

    std::shared_ptr<variantData> truth_ptr(new variantData());
    std::shared_ptr<variantData> truth_ptr2(new variantData()); // second round eval
    parse_variants(g.truth_vcf_fn, truth_ptr, truth_ptr2, ref_ptr, TRUTH);
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
    std::shared_ptr<superclusterData> sc_data_ptr(
            new superclusterData(query_ptr, truth_ptr, ref_ptr));

    // calculate supercluster sizes
    auto sc_groups = sort_superclusters(sc_data_ptr);
    g.timers[TIME_SUPCLUST].stop();

    // calculate precision/recall and genotypes
    g.timers[TIME_PR_ALN].start();
    precision_recall_threads_wrapper(sc_data_ptr, sc_groups);
    INFO("    done with precision-recall (round #1)");
    g.timers[TIME_PR_ALN].stop();

    // pull out FP and FN variants from first-round evaluation
    std::shared_ptr<variantData> query_ptr_fp(new variantData());
    std::shared_ptr<variantData> truth_ptr_fn(new variantData());
    extract_errors(sc_data_ptr, query_ptr_fp, truth_ptr_fn);

    // merge first-round FP and FN with large second-round variants
    query_ptr2->merge(query_ptr_fp);
    truth_ptr2->merge(truth_ptr_fn);

    // perform second-round evaluation
    g.timers[TIME_CLUST].start();
    simple_cluster(query_ptr2, QUERY);
    simple_cluster(truth_ptr2, TRUTH);
    g.timers[TIME_CLUST].stop();
    g.timers[TIME_SUPCLUST].start();
    std::shared_ptr<superclusterData> sc_data_ptr2(
            new superclusterData(query_ptr2, truth_ptr2, ref_ptr));
    sc_groups = sort_superclusters(sc_data_ptr2);
    g.timers[TIME_SUPCLUST].stop();
    g.timers[TIME_PR_ALN].start();
    precision_recall_threads_wrapper(sc_data_ptr2, sc_groups);
    INFO("    done with precision-recall (round #2)");
    g.timers[TIME_PR_ALN].stop();

    // merge first-round and second-round results
    // TODO: re-use add_variant_data to merge superclusters?

    // calculate global phasings
    g.timers[TIME_PHASE].start();
    std::unique_ptr<phaseblockData> phasedata_ptr(new phaseblockData(sc_data_ptr));
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
