/**
 * @file phase.h
 * @brief Phase block detection, switch/flip error classification, and phasing summary output.
 */
#ifndef _PHASE_H_
#define _PHASE_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "cluster.h"

/**
 * @class ctgPhaseblocks
 * @brief Store phase block data for one contig.
 */
class ctgPhaseblocks {
public:
    ctgPhaseblocks() {;}

    std::shared_ptr<ctgSuperclusters> ctg_superclusters = nullptr; ///< Pointer to the contig's supercluster data

    int n = 0; ///< Number of phase blocks

    // data, filled out during phase()
    std::vector<int> phase_blocks; ///< n+1 query variant indices of new phase set boundaries
    std::vector<int> switches;     ///< Variant indices before which a phase switch occurs
    std::vector<int> flips;        ///< Indices of flipped (incorrectly phased) variants
    int nswitches = 0;             ///< Number of phase switch errors
    int nflips = 0;                ///< Number of phase flip errors
};

/**
 * @class phaseblockData
 * @brief Store phase block data across all contigs.
 */
class phaseblockData {
public:
    /** @brief Constructs phaseblock container from supercluster data and runs phasing pipeline. */
    phaseblockData(std::shared_ptr<superclusterData> clusterdata_ptr);

    /**
     * @brief Writes a summary VCF containing all variants annotated with benchmark metrics.
     * @note FORMAT fields include: TP/FP/FN decision, credit score, edit distances, phase info, and flip/switch errors
     */
    void write_summary_vcf(std::string out_vcf_fn);

    /** @brief Writes phasing summary statistics to TSV file. */
    void write_phasing_summary(int phase_blocks, int switch_errors,
        int flip_errors, int variants, int ng50, int s_ngc50, int sf_ngc50);

    /**
     * @brief Writes detected switch and flip error locations and classifications to TSV file.
     * @note Error types: SWITCH_ERR, FLIP, SWITCH_AND_FLIP
     */
    void write_switchflips();

    /**
     * @brief Uses dynamic programming to find optimal phasing and detect switch/flip errors per contig.
     * @note Results stored in qvars->pb_phases and error lists in ctgPhaseblocks. Switches at phase set boundaries incur no cost.
     */
    void phase();

    /**
     * @brief Propagates phase set tags to unphased and homozygous variants.
     * @note Must run before the phase block scan, phase(), and fix_allele_counts().
     * @todo Only set phase sets for 1|1 variants when unphased evaluation is added.
     */
    void fix_phase_set_tags();

    /**
     * @brief Corrects calculated genotypes to preserve allele counts matching original calls.
     * @note Records each variant's allele count error type on both callsets, and tracks and reports
     *       genotype error statistics (0/0->0/1, 1/1->0/1, etc.)
     */
    void fix_allele_counts();

    /** @brief Writes allele count error cross-tabulation table to TSV file. */
    void write_genotype_error_summary(
            const EnumArray<ac_errtype_t, std::vector<int>, AC_ERRTYPE_SLOTS> &
                    allele_error_counts);

    /**
     * @brief Calculates NGC50 of phase blocks, optionally broken at switch or flip errors.
     * @return NGC50 value (block length at 50% cumulative length), or 0 if no blocks exist
     */
    int calculate_ng50(bool break_on_switch = false, bool break_on_flip = false);

    std::shared_ptr<fastaData> ref;   ///< Pointer to reference FASTA data
    std::vector<std::string> contigs; ///< List of all contig names
    std::vector<int> lengths;         ///< List of all contig lengths
    std::unordered_map<std::string,   ///< Map from contig name to ctgPhaseblocks
        std::shared_ptr<ctgPhaseblocks> > phase_blocks;
};

/**
 * @brief Returns the reference span of each correctly-phased block on one contig.
 * @note Phase set boundaries always break a block; switch and flip errors break one only when
 *   requested. A flip breaks twice, excising the flipped variant into a block of its own.
 */
std::vector<int> correct_block_sizes(const std::shared_ptr<ctgPhaseblocks> & ctg_pbs,
        const std::shared_ptr<ctgVariants> & qvars, bool break_on_switch, bool break_on_flip);

#endif
