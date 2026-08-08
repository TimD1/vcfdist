/**
 * @file main.cpp
 * @brief Vcfdist main entry point.
 */
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

/**
 * @brief Orchestrates the vcfdist variant benchmarking pipeline in five stages:
 *   (1) read/parse input VCFs and reference,
 *   (2) cluster variants into superclusters,
 *   (3) align and evaluate precision/recall per supercluster,
 *   (4) compute phasing statistics and detect switch/flip errors,
 *   (5) write all output files and print summary.
 * @param[in] argc Argument count
 * @param[in] argv Argument vector
 * @return EXIT_SUCCESS
 */
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);
    g.init_timers();

    g.stage(TIME_WRITE).start();
    write_params();
    g.stage(TIME_WRITE).stop();

    // parse reference fasta
    g.stage(TIME_TOTAL).start();
    g.stage(TIME_READ).start();
    std::shared_ptr<fastaData> ref_ptr(new fastaData(g.ref_fasta_fp));

    // warn before any evaluation if the stratification regions name contigs the reference lacks
    check_strata_contigs(ref_ptr);

    // parse query and truth VCFs
    std::shared_ptr<variantData> query_ptr(new variantData());
    parse_variants(g.query_vcf_fn, query_ptr, ref_ptr, QUERY);

    std::shared_ptr<variantData> truth_ptr(new variantData());
    parse_variants(g.truth_vcf_fn, truth_ptr, ref_ptr, TRUTH);
    g.stage(TIME_READ).stop();

    // ensure each input contains all contigs in BED
    intersect_contigs(query_ptr, truth_ptr, ref_ptr);

    // cluster query VCF
    g.stage(TIME_CLUSTER).start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[Q %d/%d] Clustering %s VCF%s '%s'", 
            COLOR_PURPLE, int(idx(TIME_CLUSTER)), int(idx(TIME_TOTAL))-1, callset_strs[QUERY].data(), 
            COLOR_WHITE, query_ptr->filename.data());
    std::vector<std::thread> threads;
    for (int t = 0; t < HAPS*int(query_ptr->contigs.size()); t++) {
        threads.push_back(std::thread( wf_swg_cluster, 
                    query_ptr.get(), t/2 /* contig */, static_cast<hap_t>(t%2),
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
            COLOR_PURPLE, int(idx(TIME_CLUSTER)), int(idx(TIME_TOTAL))-1, callset_strs[TRUTH].data(), 
            COLOR_WHITE, truth_ptr->filename.data());
    threads.clear();
    for (int t = 0; t < HAPS*int(truth_ptr->contigs.size()); t++) {
        threads.push_back(std::thread( wf_swg_cluster, 
                    truth_ptr.get(), t/2 /* contig */, static_cast<hap_t>(t%2),
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
            COLOR_PURPLE, int(idx(TIME_CLUSTER)), int(idx(TIME_TOTAL))-1, COLOR_WHITE);
    std::shared_ptr<superclusterData> sc_data_ptr(
            new superclusterData(query_ptr, truth_ptr, ref_ptr));
    auto sc_groups = sort_superclusters(sc_data_ptr);
    g.stage(TIME_CLUSTER).stop();

    // evaluation: precision/recall and genotypes
    g.stage(TIME_ALIGN_EVAL).start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%d/%d] Evaluating variant calls %s",
            COLOR_PURPLE, int(idx(TIME_ALIGN_EVAL)), int(idx(TIME_TOTAL))-1, COLOR_WHITE);
    precision_recall_threads_wrapper(sc_data_ptr, sc_groups);
    INFO("    done with variant call evaluation");
    g.stage(TIME_ALIGN_EVAL).stop();

    // calculate phasing statistics
    g.stage(TIME_PHASE).start();
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%d/%d] Phasing superclusters%s",
            COLOR_PURPLE, int(idx(TIME_PHASE)), int(idx(TIME_TOTAL))-1, COLOR_WHITE);
    std::unique_ptr<phaseblockData> phasedata_ptr(new phaseblockData(sc_data_ptr));
    g.stage(TIME_PHASE).stop();

    // write phasing results
    g.stage(TIME_WRITE).start();
    phasedata_ptr->write_switchflips();
    write_results(phasedata_ptr);
    phasedata_ptr->write_summary_vcf(g.out_prefix + "summary.vcf");
    g.stage(TIME_WRITE).stop();

    // report timing results
    g.stage(TIME_TOTAL).stop();
    write_runtime();
    if (g.verbosity >= 1) {
        INFO(" ")
        INFO("Timers:")
        for (stage_t t : EnumRange<stage_t, STAGE_SLOTS>{}) { g.stage(t).print(idx(t)); }
    }
    return EXIT_SUCCESS;
}
