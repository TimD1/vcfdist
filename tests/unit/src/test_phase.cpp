/**
 * @file test_phase.cpp
 * @brief Unit tests for phase.cpp: phase block statistics and switch/flip classification.
 */
#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/defs.h"
#include "../../../src/globals.h"
#include "../../../src/phase.h"
#include "test_helpers.h"

namespace {

/* Local helpers **********************************************************************************/

// Six 1-base substitutions 100bp apart, so the contig's variants span [100, 601).
std::shared_ptr<ctgVariants> six_vars() {
    std::vector<var_desc> vars;
    for (int i = 0; i < 6; i++) {
        var_desc var;
        var.pos = 100 * (i + 1);
        var.rlen = 1;
        var.ref = "A";
        var.alt = "C";
        vars.push_back(var);
    }
    return make_ctgVariants("chr1", vars);
}

// Builds phase blocks from the index vectors alone; ctg_superclusters stays null, since
// correct_block_sizes() takes the query variants directly rather than reaching through it.
std::shared_ptr<ctgPhaseblocks> make_pbs(const std::vector<int> & phase_blocks,
        const std::vector<int> & switches = {}, const std::vector<int> & flips = {}) {
    std::shared_ptr<ctgPhaseblocks> pbs(new ctgPhaseblocks());
    pbs->phase_blocks = phase_blocks;
    pbs->n = int(phase_blocks.size()) - 1;
    pbs->switches = switches;
    pbs->nswitches = int(switches.size());
    pbs->flips = flips;
    pbs->nflips = int(flips.size());
    return pbs;
}

const std::string CTG = "chr1";     ///< Contig every single-contig test builds on
const int CTG_LENGTH = 604;         ///< Contig length chosen so NG50 thresholds land on a block
const int SPACING = 100;            ///< Reference distance between consecutive built variants

/**
 * @struct ctg_input
 * @brief One contig's variants plus the metadata phaseblockData copies out of superclusterData.
 */
struct ctg_input {
    std::string ctg = CTG;                    ///< Contig name
    std::shared_ptr<ctgVariants> qvars;       ///< Query variants, or nullptr for an empty callset
    std::shared_ptr<ctgVariants> tvars;       ///< Truth variants, or nullptr for an empty callset
    int length = CTG_LENGTH;                  ///< Contig length
};

/**
 * @struct pipeline_result
 * @brief Everything one phaseblockData construction produced: the object and its INFO output.
 */
struct pipeline_result {
    std::shared_ptr<phaseblockData> data; ///< Constructed phase block data
    std::string log;                      ///< All INFO/WARN output from the construction
};

/**
 * @brief Builds query variants SPACING bases apart whose genotypes yield the requested phasings.
 *
 * phase() derives phases[i] from the orig/calc genotype pair rather than reading it, so a test
 * asks for a phasing pattern and gets the genotypes that produce it: 1|0 against 1|0 for
 * PHASE_ORIG, 1|0 against 0|1 for PHASE_SWAP, and 1|1 against 1|1 for PHASE_NONE.
 * @param[in] phases Desired phasing (PHASE_ORIG, PHASE_SWAP, or PHASE_NONE) of each variant
 * @param[in] phase_sets Phase set of each variant, or empty to place them all in phase set 1
 * @param[in] ctg Contig the variants sit on
 * @return Query variants with orig_gts, calc_gts, and phase_sets set
 */
std::shared_ptr<ctgVariants> make_qvars(const std::vector<phase_t> & phases,
        const std::vector<int> & phase_sets = {}, const std::string & ctg = CTG) {
    std::vector<var_desc> descs;
    for (size_t i = 0; i < phases.size(); i++) {
        int phase_set = phase_sets.empty() ? 1 : phase_sets[i];
        gt_t orig_gt = (phases[i] == PHASE_NONE) ? GT_ALT1_ALT1 : GT_ALT1_REF;
        descs.push_back({int(i) * SPACING, 1, TYPE_SUB, "A", "C", orig_gt, 60, phase_set});
    }
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants(ctg, descs);
    for (size_t i = 0; i < phases.size(); i++) {
        switch (phases[i]) {
            case PHASE_ORIG: qvars->calc_gts[i] = GT_ALT1_REF;  break;
            case PHASE_SWAP: qvars->calc_gts[i] = GT_REF_ALT1;  break;
            default:         qvars->calc_gts[i] = GT_ALT1_ALT1; break;
        }
    }
    return qvars;
}

/**
 * @brief Runs the phaseblockData pipeline over the given contigs, capturing its INFO output.
 *
 * phase() and fix_allele_counts() each write a TSV beneath g.out_prefix, so the prefix is pointed
 * into the caller's temporary directory; verbosity is raised so that the statistics reported only
 * through INFO are observable, and stderr is redirected so they stay out of the test log.
 * @param[in] dir Temporary directory receiving the pipeline's TSVs and captured log
 * @param[in] inputs Per-contig variants and metadata, in contig order
 * @return Constructed phase block data and its captured log
 */
pipeline_result run_pipeline(const TempDir & dir, const std::vector<ctg_input> & inputs) {
    std::vector<std::string> contigs;
    std::vector<int> lengths;
    std::vector< std::shared_ptr<ctgSuperclusters> > superclusters;
    for (const ctg_input & input : inputs) {
        contigs.push_back(input.ctg);
        lengths.push_back(input.length);
        std::shared_ptr<ctgVariants> qvars = input.qvars ? input.qvars :
                make_ctgVariants(input.ctg, {});
        std::shared_ptr<ctgVariants> tvars = input.tvars ? input.tvars :
                make_ctgVariants(input.ctg, {});
        superclusters.push_back(make_ctgSuperclusters(qvars, tvars));
    }

    g.out_prefix = dir.path() + "/";
    g.verbosity = 1;
    std::string log_fn = dir.path("pipeline.log");
    pipeline_result result;
    {
        StderrToFile redirect(log_fn);
        result.data = std::shared_ptr<phaseblockData>(new phaseblockData(
                make_superclusterData(contigs, lengths, superclusters)));
    }
    result.log = read_text(log_fn);
    return result;
}

/** @brief Runs the phaseblockData pipeline over a single contig. */
pipeline_result run_pipeline(const TempDir & dir, std::shared_ptr<ctgVariants> qvars,
        std::shared_ptr<ctgVariants> tvars = nullptr, int length = CTG_LENGTH) {
    ctg_input input;
    input.qvars = qvars;
    input.tvars = tvars;
    input.length = length;
    return run_pipeline(dir, {input});
}

/**
 * @brief Builds query variants whose middle variant carries the given genotypes and credits.
 *
 * fix_allele_counts() tie-breaks on the enclosing phase block's phasing, so the middle variant is
 * flanked by two variants that pin the block to block_phase. Every genotype pair passed here
 * classifies as PHASE_NONE, leaving the flanks in sole control of the DP, except a heterozygous
 * call that agrees with itself: that pair is PHASE_ORIG, and is only used in an unswapped block.
 * @param[in] orig_gt Original genotype of the middle variant
 * @param[in] calc_gt Calculated genotype of the middle variant
 * @param[in] block_phase Phasing (PHASE_ORIG or PHASE_SWAP) to force on the enclosing block
 * @param[in] hap1_credit Credit of the middle variant on HAP1
 * @param[in] hap2_credit Credit of the middle variant on HAP2
 * @return Query variants whose middle variant is at index 1
 */
std::shared_ptr<ctgVariants> make_ac_qvars(gt_t orig_gt, gt_t calc_gt, phase_t block_phase,
        float hap1_credit = 0, float hap2_credit = 0) {
    std::shared_ptr<ctgVariants> qvars =
            make_qvars({block_phase, PHASE_NONE, block_phase});
    qvars->orig_gts[1] = orig_gt;
    qvars->calc_gts[1] = calc_gt;
    set_hap_data(qvars, HAP1, 1, ERRTYPE_TP, 0, 60, 4, 0, hap1_credit);
    set_hap_data(qvars, HAP2, 1, ERRTYPE_TP, 0, 60, 4, 0, hap2_credit);
    return qvars;
}

/** @brief Returns the phase block data of the single contig every helper above builds. */
std::shared_ptr<ctgPhaseblocks> pbs_of(const pipeline_result & result) {
    return result.data->phase_blocks[CTG];
}

/** @brief Returns the query variants of the single contig every helper above builds. */
std::shared_ptr<ctgVariants> qvars_of(const pipeline_result & result) {
    return result.data->phase_blocks[CTG]->ctg_superclusters->callset_vars[QUERY];
}

/** @brief Returns the truth variants of the single contig every helper above builds. */
std::shared_ptr<ctgVariants> tvars_of(const pipeline_result & result) {
    return result.data->phase_blocks[CTG]->ctg_superclusters->callset_vars[TRUTH];
}

/**
 * @brief Builds truth variants SPACING bases apart with the given genotypes.
 *
 * A truth variant's calc_gt is the query genotype recovered for it by alignment, which calc_gts
 * carries in place of the reference call it is initialized to.
 * @param[in] gts Original genotype and recovered query genotype of each variant, in order
 * @return Truth variants with orig_gts and calc_gts set
 */
std::shared_ptr<ctgVariants> make_tvars(const std::vector< std::pair<uint8_t, uint8_t> > & gts) {
    std::vector<var_desc> descs;
    for (size_t i = 0; i < gts.size(); i++)
        descs.push_back({int(i) * SPACING, 1, TYPE_SUB, "A", "C", gts[i].first, 60, 1});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG, descs);
    for (size_t i = 0; i < gts.size(); i++) tvars->calc_gts[i] = gts[i].second;
    return tvars;
}

/* correct_block_sizes ****************************************************************************/

TEST(CorrectBlockSizes, NoVariants) {
    // an unphased contig contributes no blocks at all, not a single zero-length one
    auto qvars = make_ctgVariants("chr1", {});
    auto pbs = make_pbs({0});
    EXPECT_EQ(std::vector<int>(), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, SingleBlockSpansAllVariants) {
    // one phase set with no errors: the block runs from the first variant to the last variant's end
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 6});
    EXPECT_EQ(std::vector<int>({501}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, PhaseSetBoundaryAlwaysBreaks) {
    // a new phase set is not an error, so it splits blocks whether or not breaking is requested
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 3, 6});
    EXPECT_EQ(std::vector<int>({201, 201}), correct_block_sizes(pbs, qvars, false, false));
    EXPECT_EQ(std::vector<int>({201, 201}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, SwitchErrorBreaksOnlyWhenRequested) {
    // the switch before variant 2 ends one block at variant 1 and starts the next at variant 2
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 6}, {2});
    EXPECT_EQ(std::vector<int>({501}), correct_block_sizes(pbs, qvars, false, false));
    EXPECT_EQ(std::vector<int>({101, 301}), correct_block_sizes(pbs, qvars, true, false));
}

TEST(CorrectBlockSizes, FlipBreaksTwiceOnlyWhenRequested) {
    // a flipped variant is wrong in isolation, so it is excised into a block of its own
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 6}, {}, {2});
    EXPECT_EQ(std::vector<int>({501}), correct_block_sizes(pbs, qvars, true, false));
    EXPECT_EQ(std::vector<int>({101, 1, 201}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, ConsecutiveFlips) {
    // each flip is excised independently, leaving the spans between them as their own blocks
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 6}, {}, {1, 3});
    EXPECT_EQ(std::vector<int>({1, 1, 1, 1, 101}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, SwitchAndFlipAdvancesBothIndices) {
    // one variant carrying both a new phase set and a flip is consumed once: the phase set
    // boundary is not revisited on the next iteration, so this yields three blocks, not four
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 2, 6}, {}, {2});
    EXPECT_EQ(std::vector<int>({101, 1, 201}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, FlipTakesPrecedenceOverPhaseSetBoundary) {
    // the flip is checked last, so a tie at the same variant is classified as switch-and-flip and
    // breaks twice; a phase set boundary alone at that variant would break only once
    auto qvars = six_vars();
    auto tie = make_pbs({0, 2, 6}, {}, {2});
    auto boundary_only = make_pbs({0, 2, 6});
    EXPECT_EQ(3u, correct_block_sizes(tie, qvars, true, true).size());
    EXPECT_EQ(2u, correct_block_sizes(boundary_only, qvars, true, true).size());
}

TEST(CorrectBlockSizes, FinalBlockAppendedAfterLastBreak) {
    // the trailing block is appended after the loop, so the last variant is never dropped
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 5, 6});
    EXPECT_EQ(std::vector<int>({401, 1}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, FlipOnLastVariant) {
    // the last variant is excised into its own block, and nothing is left to open a block after it
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 6}, {}, {5});
    EXPECT_EQ(std::vector<int>({401, 1}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, FlipOnLastVariantWithPendingSwitch) {
    // the flip takes precedence at a shared index, and the switch error left pending behind it is
    // not revisited: no block is reopened after the last variant, whatever else points at it
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 6}, {5}, {5});
    EXPECT_EQ(std::vector<int>({401, 1}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, VariableLengthVariants) {
    // block bounds come from poss and rlens, so a deletion extends the block past its start
    std::vector<var_desc> vars;
    var_desc del;
    del.pos = 100;
    del.rlen = 50;
    del.type = TYPE_DEL;
    del.ref = std::string(50, 'A');
    del.alt = "A";
    vars.push_back(del);
    var_desc sub;
    sub.pos = 300;
    sub.rlen = 1;
    sub.ref = "A";
    sub.alt = "C";
    vars.push_back(sub);
    auto qvars = make_ctgVariants("chr1", vars);
    auto pbs = make_pbs({0, 2});
    EXPECT_EQ(std::vector<int>({201}), correct_block_sizes(pbs, qvars, true, true));
}

TEST(CorrectBlockSizes, OutOfOrderIndicesError) {
    // breaks are consumed in ascending order, so a descending index list cannot be walked
    auto qvars = six_vars();
    auto pbs = make_pbs({0, 6}, {}, {3, 1});
    EXPECT_EXIT(correct_block_sizes(pbs, qvars, true, true), testing::ExitedWithCode(1),
            "is not after current variant");
}

/* phase(): variant classification ****************************************************************/

TEST(Phase, ClassifyOrig) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG}));
    EXPECT_EQ(PHASE_ORIG, qvars_of(result)->phases[0]);
}

TEST(Phase, ClassifySwap) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_SWAP}));
    EXPECT_EQ(PHASE_SWAP, qvars_of(result)->phases[0]);

    // a lone swapped variant is a self-consistent block, not an error
    EXPECT_EQ(PHASE_SWAP, qvars_of(result)->pb_phases[0]);
    EXPECT_EQ(0, pbs_of(result)->nswitches);
    EXPECT_EQ(0, pbs_of(result)->nflips);
}

TEST(Phase, ClassifyNone) {
    GlobalsGuard guard;
    TempDir dir;

    // a homozygous variant carries no phasing information either way
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_NONE}));
    EXPECT_EQ(PHASE_NONE, qvars_of(result)->phases[0]);
}

/* phase(): switch and flip errors ****************************************************************/

TEST(Phase, PerfectBlock) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}));
    EXPECT_EQ(0, pbs_of(result)->nswitches);
    EXPECT_EQ(0, pbs_of(result)->nflips);
    EXPECT_EQ(std::vector<phase_t>({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}), qvars_of(result)->pb_phases);
}

TEST(Phase, AllSwapBlock) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_SWAP, PHASE_SWAP, PHASE_SWAP}));

    // a uniformly swapped block is correctly phased relative to itself: the backtrace starts in
    // PHASE_SWAP and stays there, so no variant is reported as an error
    EXPECT_EQ(0, pbs_of(result)->nswitches);
    EXPECT_EQ(0, pbs_of(result)->nflips);
    EXPECT_EQ(std::vector<phase_t>({PHASE_SWAP, PHASE_SWAP, PHASE_SWAP}), qvars_of(result)->pb_phases);
}

TEST(Phase, SingleFlip) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_SWAP, PHASE_ORIG}));
    EXPECT_EQ(1, pbs_of(result)->nflips);
    EXPECT_EQ(std::vector<int>({1}), pbs_of(result)->flips);
    EXPECT_EQ(0, pbs_of(result)->nswitches);
}

TEST(Phase, SingleSwitch) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP, PHASE_SWAP}));

    // switches are recorded before the variant the phase changes at
    EXPECT_EQ(1, pbs_of(result)->nswitches);
    EXPECT_EQ(std::vector<int>({2}), pbs_of(result)->switches);
    EXPECT_EQ(0, pbs_of(result)->nflips);
}

TEST(Phase, SwitchPlusFlip) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG,
            PHASE_SWAP, PHASE_ORIG, PHASE_SWAP, PHASE_SWAP}));

    // one switch and one flip, tallied independently; the DP is indifferent between switching at
    // index 3 with a flip at index 4 and the reverse, and the backtrace breaks the tie this way
    EXPECT_EQ(1, pbs_of(result)->nswitches);
    EXPECT_EQ(std::vector<int>({5}), pbs_of(result)->switches);
    EXPECT_EQ(1, pbs_of(result)->nflips);
    EXPECT_EQ(std::vector<int>({3}), pbs_of(result)->flips);
}

TEST(Phase, BoundaryFree) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP, PHASE_SWAP}, {1, 1, 2, 2}));

    // the phasing changes across a phase set boundary, where a swap costs nothing and is not an
    // error: the two phase sets are each internally consistent
    EXPECT_EQ(0, pbs_of(result)->nswitches);
    EXPECT_EQ(0, pbs_of(result)->nflips);
    EXPECT_EQ(2, pbs_of(result)->n);
    EXPECT_EQ(std::vector<phase_t>({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP, PHASE_SWAP}),
            qvars_of(result)->pb_phases);
}

TEST(Phase, MultipleBlocks) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_SWAP, PHASE_ORIG, PHASE_SWAP, PHASE_ORIG, PHASE_SWAP},
                       {1, 1, 1, 2, 2, 2}));

    // each phase set resolves to its own majority phasing, leaving one flip in each
    EXPECT_EQ(0, pbs_of(result)->nswitches);
    EXPECT_EQ(2, pbs_of(result)->nflips);
    EXPECT_EQ(std::vector<int>({1, 4}), pbs_of(result)->flips);
    EXPECT_EQ(std::vector<phase_t>({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG,
            PHASE_SWAP, PHASE_SWAP, PHASE_SWAP}), qvars_of(result)->pb_phases);
}

TEST(Phase, EmptyContig) {
    GlobalsGuard guard;
    TempDir dir;

    // the backwards pass is skipped entirely, leaving the tallies at zero
    pipeline_result result = run_pipeline(dir, nullptr);
    EXPECT_EQ(0, pbs_of(result)->nswitches);
    EXPECT_EQ(0, pbs_of(result)->nflips);
    EXPECT_TRUE(pbs_of(result)->switches.empty());
    EXPECT_TRUE(pbs_of(result)->flips.empty());
}

TEST(Phase, NoneNotFlipped) {
    GlobalsGuard guard;
    TempDir dir;

    // an unphased variant between two agreeing variants costs nothing in either phasing, so it
    // cannot disagree with the block and is never reported as a flip
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_NONE, PHASE_ORIG}));
    EXPECT_EQ(0, pbs_of(result)->nflips);
    EXPECT_EQ(0, pbs_of(result)->nswitches);
}

TEST(Phase, FlipsSwitchesAscending) {
    GlobalsGuard guard;
    TempDir dir;

    // three runs of five, each with a flipped variant in the middle
    pipeline_result result = run_pipeline(dir, make_qvars({
            PHASE_ORIG, PHASE_ORIG, PHASE_SWAP, PHASE_ORIG, PHASE_ORIG,
            PHASE_SWAP, PHASE_SWAP, PHASE_ORIG, PHASE_SWAP, PHASE_SWAP,
            PHASE_ORIG, PHASE_ORIG, PHASE_SWAP, PHASE_ORIG, PHASE_ORIG}));

    // both lists are built back-to-front and reversed, so they read in ascending index order
    EXPECT_EQ(std::vector<int>({5, 10}), pbs_of(result)->switches);
    EXPECT_EQ(std::vector<int>({2, 7, 12}), pbs_of(result)->flips);
    EXPECT_TRUE(std::is_sorted(pbs_of(result)->switches.begin(), pbs_of(result)->switches.end()));
    EXPECT_TRUE(std::is_sorted(pbs_of(result)->flips.begin(), pbs_of(result)->flips.end()));
    EXPECT_EQ(size_t(pbs_of(result)->nswitches), pbs_of(result)->switches.size());
    EXPECT_EQ(size_t(pbs_of(result)->nflips), pbs_of(result)->flips.size());
}

TEST(Phase, PbPhasesAssignment) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP, PHASE_SWAP}));

    // the backtrace writes pb_phases[i-1] after following the pointer at i, so the switch recorded
    // at index 2 is the first index carrying the swapped block phasing
    EXPECT_EQ(size_t(4), qvars_of(result)->pb_phases.size());
    EXPECT_EQ(std::vector<phase_t>({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP, PHASE_SWAP}),
            qvars_of(result)->pb_phases);
}

/* fix_phase_set_tags() ***************************************************************************/

TEST(FixPhaseSetTags, NoPhaseSetsContig) {
    GlobalsGuard guard;
    TempDir dir;

    // no variant is phased, so there is nothing to propagate and the tags are left alone
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG}, {0, 0}));
    EXPECT_EQ(std::vector<int>({0, 0}), qvars_of(result)->phase_sets);
}

TEST(FixPhaseSetTags, BackfillLeadingZeros) {
    GlobalsGuard guard;
    TempDir dir;

    // leading unphased variants join the first phase set that appears
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {0, 0, 5, 5}));
    EXPECT_EQ(std::vector<int>({5, 5, 5, 5}), qvars_of(result)->phase_sets);
}

TEST(FixPhaseSetTags, PropagateUnphasedMiddle) {
    GlobalsGuard guard;
    TempDir dir;

    // an unphased variant joins the phase set of the variants before it
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {7, 0, 7}));
    EXPECT_EQ(std::vector<int>({7, 7, 7}), qvars_of(result)->phase_sets);
}

TEST(FixPhaseSetTags, NewPsSpanReset) {
    GlobalsGuard guard;
    TempDir dir;

    // two phase sets of two variants each: spans are 0-101 and 200-301, both 101 bases wide, and
    // the second span starts over rather than continuing from the first. Spans reach the caller
    // only through the reported NG50, which is 101 for {101, 101} against 202 total bases; a span
    // continuing from the first phase set would give {101, 301} and an NG50 of 301.
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {1, 1, 2, 2}), nullptr,
            202);
    EXPECT_TRUE(logged(result.log, "QUERY phase block NG50: 101")) << result.log;
}

TEST(FixPhaseSetTags, SamePsExtendsEnd) {
    GlobalsGuard guard;
    TempDir dir;

    // the middle variant is a 201-base deletion reaching past the end of the last variant, so the
    // running span end is a maximum rather than the last variant's end
    std::shared_ptr<ctgVariants> qvars = make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG});
    qvars->rlens[1] = 201;
    pipeline_result result = run_pipeline(dir, qvars, nullptr, 602);

    // the span is 0-301, so NG50 is 301 against 602 total bases; taking the last variant's end
    // instead of the running maximum would give 201
    EXPECT_TRUE(logged(result.log, "QUERY phase block NG50: 301")) << result.log;
}

TEST(FixPhaseSetTags, FinalSpanPushed) {
    GlobalsGuard guard;
    TempDir dir;

    // a lone phase set is never followed by another, so the push inside the variant loop never
    // fires: the span 0-201 is recorded only by the push after the loop, and without it there
    // would be no block at all and the reported NG50 would be 0
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}), nullptr, 402);
    EXPECT_TRUE(logged(result.log, "QUERY phase block NG50: 201")) << result.log;
}

TEST(FixPhaseSetTags, PhaseSetCountMatchesPhaseSets) {
    GlobalsGuard guard;
    TempDir dir;

    // two phase sets are reported as two; the empty truth callset takes the early exit, which
    // contributes the whole contig as one span and so reports one
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {1, 1, 2, 2}));
    EXPECT_TRUE(logged(result.log, "QUERY phase sets: 2")) << result.log;
    EXPECT_TRUE(logged(result.log, "TRUTH phase sets: 1")) << result.log;
}

TEST(FixPhaseSetTags, PhaseSetCountSumsAcrossContigs) {
    GlobalsGuard guard;
    TempDir dir;

    // the count accumulates across contigs, so an error made once per phased contig compounds:
    // two contigs of two phase sets each are reported as four, not six
    ctg_input first;
    first.qvars = make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {1, 1, 2, 2});
    ctg_input second;
    second.ctg = "chr2";
    second.qvars = make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {3, 3, 4, 4},
            "chr2");
    pipeline_result result = run_pipeline(dir, {first, second});
    EXPECT_TRUE(logged(result.log, "QUERY phase sets: 4")) << result.log;
}

TEST(FixPhaseSetTags, BothCallsets) {
    GlobalsGuard guard;
    TempDir dir;

    // truth tags are propagated on the same pass as query tags
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG,
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 0},
             {100, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 0},
             {200, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 3}});
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG}, {0, 9}), tvars);
    EXPECT_EQ(std::vector<int>({9, 9}), qvars_of(result)->phase_sets);
    EXPECT_EQ(std::vector<int>({3, 3, 3}),
            pbs_of(result)->ctg_superclusters->callset_vars[TRUTH]->phase_sets);
}

TEST(FixPhaseSetTags, Ng50Reported) {
    GlobalsGuard guard;
    TempDir dir;

    // phase set count, NG50, and total bases are reported for both callsets
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG, PHASE_ORIG}));
    EXPECT_TRUE(logged(result.log, "QUERY phase sets:"));
    EXPECT_TRUE(logged(result.log, "QUERY phase block NG50:"));
    EXPECT_TRUE(logged(result.log, "QUERY total bases:"));
    EXPECT_TRUE(logged(result.log, "TRUTH phase sets:"));
    EXPECT_TRUE(logged(result.log, "TRUTH phase block NG50:"));
    EXPECT_TRUE(logged(result.log, "TRUTH total bases:"));
}

TEST(FixPhaseSetTags, UnphasedContigsEachContributeOwnLength) {
    GlobalsGuard guard;
    TempDir dir;

    // each contig with no phase sets contributes its own whole length as one span, so the spans are
    // {100, 500} against 600 total bases and the NG50 is 500. Reusing the first contig's length for
    // the second would give {100, 100}, which never reaches half the genome and reports 0.
    ctg_input first;
    first.ctg = "chr1";
    first.qvars = make_qvars({PHASE_ORIG, PHASE_ORIG}, {0, 0});
    first.length = 100;
    ctg_input second;
    second.ctg = "chr2";
    second.qvars = make_ctgVariants("chr2", {});
    second.length = 500;
    pipeline_result result = run_pipeline(dir, {first, second});
    EXPECT_TRUE(logged(result.log, "QUERY phase block NG50: 500")) << result.log;
}

TEST(FixPhaseSetTags, UnphasedContigLengthNotTakenFromPhasedContig) {
    GlobalsGuard guard;
    TempDir dir;

    // an unphased contig looks up its own length even with a phased contig in between: the spans
    // are chr1's 100 bases, chr2's phased 0-101, and chr3's 1000 bases, so the NG50 is 1000 against
    // 1400 total bases. Counting only the phased contigs would hand chr3 chr2's 300 bases, and the
    // resulting {100, 101, 300} never reaches half the genome and reports 0.
    ctg_input first;
    first.ctg = "chr1";
    first.qvars = make_qvars({PHASE_ORIG, PHASE_ORIG}, {0, 0});
    first.length = 100;
    ctg_input second;
    second.ctg = "chr2";
    second.qvars = make_qvars({PHASE_ORIG, PHASE_ORIG}, {}, "chr2");
    second.length = 300;
    ctg_input third;
    third.ctg = "chr3";
    third.qvars = make_ctgVariants("chr3", {});
    third.length = 1000;
    pipeline_result result = run_pipeline(dir, {first, second, third});
    EXPECT_TRUE(logged(result.log, "QUERY phase block NG50: 1000")) << result.log;
}

/* fix_allele_counts() ****************************************************************************/

TEST(FixAlleleCounts, UnknownErrors) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";

    // a reference call that stayed a reference call has no allele count error type
    std::shared_ptr<ctgVariants> qvars = make_qvars({PHASE_NONE});
    qvars->orig_gts[0] = GT_REF_REF;
    qvars->calc_gts[0] = GT_REF_REF;
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {CTG}, {CTG_LENGTH}, {make_ctgSuperclusters(qvars, make_ctgVariants(CTG, {}))});
    EXPECT_EXIT(phaseblockData data(sc_data), testing::ExitedWithCode(1),
            "Unknown variant allele count");
}

TEST(FixAlleleCounts, OneToOneTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_REF_ALT1, PHASE_ORIG));
    EXPECT_EQ(AC_ERR_1_TO_1, qvars_of(result)->ac_errtype[1]);
    EXPECT_EQ(GT_REF_ALT1, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, OneToTwoTallied) {
    GlobalsGuard guard;
    TempDir dir;

    // a heterozygous truth allele called homozygous: one allele is right, the other spurious
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_ALT1_ALT1, GT_REF_ALT1, PHASE_ORIG));
    EXPECT_EQ(AC_ERR_1_TO_2, qvars_of(result)->ac_errtype[1]);
}

TEST(FixAlleleCounts, ZeroToTwoTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_ALT1_ALT1, GT_REF_REF, PHASE_ORIG));
    EXPECT_EQ(AC_ERR_0_TO_2, qvars_of(result)->ac_errtype[1]);
}

// The truth callset gets the same value the query callset does, read from the other side: a truth
// record's own orig_gt supplies the truth allele count and its recovered calc_gt the query's.

TEST(FixAlleleCounts, TruthHeterozygousFalseNegativeTallied) {
    GlobalsGuard guard;
    TempDir dir;

    // one truth allele that no query allele matched, so its recovered genotype stayed 0|0
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT1, GT_REF_REF}}));
    EXPECT_EQ(AC_ERR_1_TO_0, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, TruthHomozygousFalseNegativeTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_ALT1_ALT1, GT_REF_REF}}));
    EXPECT_EQ(AC_ERR_2_TO_0, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, TruthGenotypeErrorTallied) {
    GlobalsGuard guard;
    TempDir dir;

    // one truth allele the query called on both haplotypes: the truth record reports it too
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT1, GT_ALT1_ALT1}}));
    EXPECT_EQ(AC_ERR_1_TO_2, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, TruthHomozygousMatchTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_ALT1_ALT1, GT_ALT1_ALT1}}));
    EXPECT_EQ(AC_ERR_2_TO_2, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, NoTruthVariantLeftUnknown) {
    GlobalsGuard guard;
    TempDir dir;

    // a pure false negative, a genotype error, and a homozygous match all reach a defined value
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT1, GT_REF_REF},
                        {GT_REF_ALT1, GT_ALT1_ALT1},
                        {GT_ALT1_ALT1, GT_ALT1_ALT1}}));
    for (int vi = 0; vi < tvars_of(result)->n; vi++)
        EXPECT_NE(AC_UNKNOWN, tvars_of(result)->ac_errtype[vi]) << "truth variant " << vi;
}

TEST(FixAlleleCounts, TruthValuesStayOutOfTheGenotypeErrorSummary) {
    GlobalsGuard guard;
    TempDir dir;

    // the summary tallies the query loop plus the hand-rolled truth false-negative branches, so a
    // truth record now reporting 0/1 -> 1/1 must not add a second count to that row
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT1, GT_ALT1_ALT1}}));
    EXPECT_TRUE(logged(result.log, "0/1 -> 1/1: 0")) << result.log;
}

TEST(FixAlleleCounts, ForceOneOneKeepsGt) {
    GlobalsGuard guard;
    TempDir dir;

    // a 1|1 call is always evaluated as 1|1, whatever the alignment calculated
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_ALT1_ALT1, GT_ALT1_REF, PHASE_ORIG));
    EXPECT_EQ(GT_ALT1_ALT1, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, ForceOneOneSwapsHapData) {
    GlobalsGuard guard;
    TempDir dir;

    // a calc_gt of 1|1 carries no record of which haplotype its data came from, so inside a
    // swapped block the per-haplotype lanes are exchanged when the genotype is forced back
    std::shared_ptr<ctgVariants> qvars = make_ac_qvars(GT_ALT1_ALT1, GT_ALT1_ALT1, PHASE_SWAP);
    set_hap_data(qvars, HAP1, 1, ERRTYPE_TP, 1, 10, 2, 3, 0.25);
    set_hap_data(qvars, HAP2, 1, ERRTYPE_FP, 2, 20, 4, 5, 0.75);
    pipeline_result result = run_pipeline(dir, qvars);

    ASSERT_EQ(PHASE_SWAP, qvars_of(result)->pb_phases[1]);
    EXPECT_EQ(ERRTYPE_FP, qvars_of(result)->errtypes[HAP1][1]);
    EXPECT_EQ(ERRTYPE_TP, qvars_of(result)->errtypes[HAP2][1]);
    EXPECT_EQ(2, qvars_of(result)->sync_group[HAP1][1]);
    EXPECT_EQ(1, qvars_of(result)->sync_group[HAP2][1]);
    EXPECT_FLOAT_EQ(20, qvars_of(result)->callq[HAP1][1]);
    EXPECT_FLOAT_EQ(10, qvars_of(result)->callq[HAP2][1]);
    EXPECT_EQ(4, qvars_of(result)->ref_ed[HAP1][1]);
    EXPECT_EQ(2, qvars_of(result)->ref_ed[HAP2][1]);
    EXPECT_EQ(5, qvars_of(result)->query_ed[HAP1][1]);
    EXPECT_EQ(3, qvars_of(result)->query_ed[HAP2][1]);
    EXPECT_FLOAT_EQ(0.75, qvars_of(result)->credit[HAP1][1]);
    EXPECT_FLOAT_EQ(0.25, qvars_of(result)->credit[HAP2][1]);
}

TEST(FixAlleleCounts, TwoToOneHap1Better) {
    GlobalsGuard guard;
    TempDir dir;

    // a heterozygous call evaluated as homozygous keeps the haplotype with better credit
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_ALT1_ALT1, PHASE_ORIG, 0.9, 0.1));
    EXPECT_EQ(AC_ERR_2_TO_1, qvars_of(result)->ac_errtype[1]);
    EXPECT_EQ(GT_ALT1_REF, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, TwoToOneHap2Better) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_ALT1_ALT1, PHASE_ORIG, 0.1, 0.9));
    EXPECT_EQ(GT_REF_ALT1, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, TwoToOneTieOrig) {
    GlobalsGuard guard;
    TempDir dir;

    // equal credit on both haplotypes falls back to the block's phasing, which keeps the call
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_ALT1_ALT1, PHASE_ORIG, 0.5, 0.5));
    EXPECT_EQ(GT_REF_ALT1, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, TwoToOneTieSwap) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_ALT1_ALT1, PHASE_SWAP, 0.5, 0.5));
    ASSERT_EQ(PHASE_SWAP, qvars_of(result)->pb_phases[1]);
    EXPECT_EQ(GT_ALT1_REF, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneHap1Better) {
    GlobalsGuard guard;
    TempDir dir;

    // a heterozygous call evaluated as reference is restored on the better-credit haplotype
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_REF_REF, PHASE_ORIG, 0.9, 0.1));
    EXPECT_EQ(AC_ERR_0_TO_1, qvars_of(result)->ac_errtype[1]);
    EXPECT_EQ(GT_ALT1_REF, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneHap2Better) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_REF_REF, PHASE_ORIG, 0.1, 0.9));
    EXPECT_EQ(GT_REF_ALT1, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneTieOrig) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_REF_REF, PHASE_ORIG, 0.5, 0.5));
    EXPECT_EQ(GT_REF_ALT1, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneTieSwap) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT1, GT_REF_REF, PHASE_SWAP, 0.5, 0.5));
    ASSERT_EQ(PHASE_SWAP, qvars_of(result)->pb_phases[1]);
    EXPECT_EQ(GT_ALT1_REF, qvars_of(result)->calc_gts[1]);
}

TEST(FixAlleleCounts, TruthFn2To0) {
    GlobalsGuard guard;
    TempDir dir;

    // a homozygous truth variant missed on both haplotypes loses two alleles, and only the truth
    // VCF can show it: there is no query variant to tally it against
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG,
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT1_ALT1, 60, 1}});
    set_hap_data(tvars, HAP1, 0, ERRTYPE_FN, 0, 0, 1, 1, 0);
    set_hap_data(tvars, HAP2, 0, ERRTYPE_FN, 0, 0, 1, 1, 0);
    pipeline_result result = run_pipeline(dir, nullptr, tvars);
    EXPECT_TRUE(logged(result.log, "1/1 -> 0/0: 1 "));
}

TEST(FixAlleleCounts, TruthFn1To0Hap1) {
    GlobalsGuard guard;
    TempDir dir;

    // for a 1|0 truth variant only HAP1 carries the allele, so only HAP1 can be a false negative
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG,
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 1}});
    set_hap_data(tvars, HAP1, 0, ERRTYPE_FN, 0, 0, 1, 1, 0);
    pipeline_result result = run_pipeline(dir, nullptr, tvars);
    EXPECT_TRUE(logged(result.log, "0/1 -> 0/0: 1 "));
}

TEST(FixAlleleCounts, TruthFn1To0Hap2) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG,
            {{0, 1, TYPE_SUB, "A", "C", GT_REF_ALT1, 60, 1}});
    set_hap_data(tvars, HAP2, 0, ERRTYPE_FN, 0, 0, 1, 1, 0);
    pipeline_result result = run_pipeline(dir, nullptr, tvars);
    EXPECT_TRUE(logged(result.log, "0/1 -> 0/0: 1 "));
}

TEST(FixAlleleCounts, TruthLoopBound) {
    GlobalsGuard guard;
    TempDir dir;

    // the last truth variant is inside the loop bound (#67): it is the only false negative here,
    // so a bound one short would report none
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG,
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 1},
             {100, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 1},
             {200, 1, TYPE_SUB, "A", "C", GT_ALT1_ALT1, 60, 1}});
    set_hap_data(tvars, HAP1, 0, ERRTYPE_TP, 0, 60, 1, 0, 1);
    set_hap_data(tvars, HAP1, 1, ERRTYPE_TP, 0, 60, 1, 0, 1);
    set_hap_data(tvars, HAP1, 2, ERRTYPE_FN, 0, 0, 1, 1, 0);
    set_hap_data(tvars, HAP2, 2, ERRTYPE_FN, 0, 0, 1, 1, 0);
    pipeline_result result = run_pipeline(dir, nullptr, tvars);
    EXPECT_TRUE(logged(result.log, "1/1 -> 0/0: 1 "));
}

/* calculate_ng50() *******************************************************************************/

TEST(PhaseblockNg50, NoBreaks) {
    GlobalsGuard guard;
    TempDir dir;

    // six variants 100 bases apart in one phase set: one block spanning 0-501
    pipeline_result result = run_pipeline(dir, make_qvars(std::vector<phase_t>(6, PHASE_ORIG)));
    EXPECT_EQ(501, result.data->calculate_ng50(false, false));
}

TEST(PhaseblockNg50, SwitchBreak) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG,
            PHASE_SWAP, PHASE_SWAP, PHASE_SWAP}));
    ASSERT_EQ(std::vector<int>({3}), pbs_of(result)->switches);

    // the switch splits 0-501 into 0-201 and 300-501
    EXPECT_EQ(201, result.data->calculate_ng50(true, false));
}

TEST(PhaseblockNg50, SwitchflipBreak) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP,
            PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}));
    ASSERT_EQ(std::vector<int>({2}), pbs_of(result)->flips);

    // a flip breaks twice, isolating the flipped variant: 0-101, the variant itself, and 300-501
    EXPECT_EQ(101, result.data->calculate_ng50(true, true));
}

TEST(PhaseblockNg50, FlipIgnoredWhenOff) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP,
            PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}));
    ASSERT_EQ(1, pbs_of(result)->nflips);

    // the same flip leaves the block whole when only switches break it
    EXPECT_EQ(501, result.data->calculate_ng50(true, false));
}

TEST(PhaseblockNg50, SwitchIgnoredWhenOff) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG,
            PHASE_SWAP, PHASE_SWAP, PHASE_SWAP}));
    ASSERT_EQ(1, pbs_of(result)->nswitches);

    // NG50 rather than NGC50: the switch error is not a break
    EXPECT_EQ(501, result.data->calculate_ng50(false, false));
}

TEST(PhaseblockNg50, FlipOnLastVariant) {
    GlobalsGuard guard;
    TempDir dir;

    // phase() records a trailing swapped variant as a flip on the last variant, so this is the
    // shortest input reaching the last-variant flip through the pipeline rather than by hand
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_SWAP}), nullptr, 202);
    ASSERT_EQ(std::vector<int>({2}), pbs_of(result)->flips);

    // the flip leaves 0-101 and the flipped variant itself, and no block beyond it
    EXPECT_EQ(101, result.data->calculate_ng50(true, true));
}

TEST(PhaseblockNg50, EmptyReturnsZero) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, nullptr);
    EXPECT_EQ(0, result.data->calculate_ng50(false, false));
}

TEST(PhaseblockNg50, ThresholdSelection) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_qvars(std::vector<phase_t>(6, PHASE_ORIG)));

    // blocks of 201, 101, and 1 bases against 604 total: the running sum first reaches half the
    // genome at the second block, so the second block's length is the answer
    pbs_of(result)->switches = {3, 5};
    pbs_of(result)->nswitches = 2;
    EXPECT_EQ(101, result.data->calculate_ng50(true, false));
}

TEST(PhaseblockNg50, NeverReachesHalfZero) {
    GlobalsGuard guard;
    TempDir dir;

    // 501 phased bases out of 10000 never reach half the genome
    pipeline_result result = run_pipeline(dir, make_qvars(std::vector<phase_t>(6, PHASE_ORIG)),
            nullptr, 10000);
    EXPECT_EQ(0, result.data->calculate_ng50(false, false));
}

TEST(PhaseblockNg50, MultiContigBlocksPooled) {
    GlobalsGuard guard;
    TempDir dir;

    // correct_block_sizes() works one contig at a time; the reduction is over every contig's
    // blocks at once, against the whole genome. chr1's 501-base block falls one base short of
    // half of 1004, so the running sum only crosses half at chr2's 101-base block — which it can
    // only reach if both contigs' blocks are in one list. With chr1's alone the answer would be 0.
    ctg_input first;
    first.ctg = "chr1";
    first.qvars = make_qvars(std::vector<phase_t>(6, PHASE_ORIG));
    first.length = 502;
    ctg_input second;
    second.ctg = "chr2";
    second.qvars = make_ctgVariants("chr2",
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 1},
             {100, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 1}});
    second.qvars->calc_gts[0] = GT_ALT1_REF;
    second.qvars->calc_gts[1] = GT_ALT1_REF;
    second.length = 502;
    pipeline_result result = run_pipeline(dir, {first, second});
    EXPECT_EQ(101, result.data->calculate_ng50(false, false));
}

/* phaseblockData constructor *********************************************************************/

TEST(PhaseblockDataCtor, CopiesMetadata) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<fastaData> ref = make_fasta({{"chr1", std::string(300, 'A')},
            {"chr2", std::string(300, 'C')}});
    g.out_prefix = dir.path() + "/";
    g.verbosity = 0;

    std::shared_ptr<ctgSuperclusters> chr1_scs = make_ctgSuperclusters(
            make_qvars({PHASE_ORIG, PHASE_ORIG}), make_ctgVariants("chr1", {}));
    std::shared_ptr<ctgSuperclusters> chr2_scs = make_ctgSuperclusters(
            make_ctgVariants("chr2", {}), make_ctgVariants("chr2", {}));
    phaseblockData data(make_superclusterData(
            {"chr1", "chr2"}, {300, 250}, {chr1_scs, chr2_scs}, ref));

    EXPECT_EQ(std::vector<std::string>({"chr1", "chr2"}), data.contigs);
    EXPECT_EQ(std::vector<int>({300, 250}), data.lengths);
    EXPECT_EQ(ref, data.ref);
    EXPECT_EQ(size_t(2), data.phase_blocks.size());
    EXPECT_EQ(chr1_scs, data.phase_blocks["chr1"]->ctg_superclusters);
    EXPECT_EQ(chr2_scs, data.phase_blocks["chr2"]->ctg_superclusters);
}

TEST(PhaseblockDataCtor, BoundariesSinglePs) {
    GlobalsGuard guard;
    TempDir dir;

    // boundaries bracket the variants: one block, from index 0 up to the variant count
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}));
    EXPECT_EQ(1, pbs_of(result)->n);
    EXPECT_EQ(std::vector<int>({0, 3}), pbs_of(result)->phase_blocks);
}

TEST(PhaseblockDataCtor, BoundariesMultiPs) {
    GlobalsGuard guard;
    TempDir dir;

    // a boundary is recorded at every index whose phase set differs from the one before it
    pipeline_result result = run_pipeline(dir, make_qvars(
            {PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {1, 1, 2, 2, 2}));
    EXPECT_EQ(2, pbs_of(result)->n);
    EXPECT_EQ(std::vector<int>({0, 2, 5}), pbs_of(result)->phase_blocks);
}

TEST(PhaseblockDataCtor, BoundaryAtIndexZero) {
    GlobalsGuard guard;
    TempDir dir;

    // the scan starts from a phase set of -1, so the first variant always opens a block, even
    // when its phase set tag is missing
    pipeline_result result = run_pipeline(dir, make_qvars({PHASE_ORIG, PHASE_ORIG}, {0, 0}));
    EXPECT_EQ(1, pbs_of(result)->n);
    EXPECT_EQ(std::vector<int>({0, 2}), pbs_of(result)->phase_blocks);
}

TEST(PhaseblockDataCtor, EmptyQueryContig) {
    GlobalsGuard guard;
    TempDir dir;

    // no variants, so only the trailing boundary is pushed and no block is counted
    pipeline_result result = run_pipeline(dir, nullptr);
    EXPECT_EQ(0, pbs_of(result)->n);
    EXPECT_EQ(std::vector<int>({0}), pbs_of(result)->phase_blocks);
}

TEST(PhaseblockDataCtor, PipelineOrder) {
    GlobalsGuard guard;
    TempDir dir;

    // boundaries are scanned after fix_phase_set_tags() backfills the missing tags, so the leading
    // unphased pair joins the phase set it is backfilled into rather than opening a block of its own
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {0, 0, 5, 5}));
    EXPECT_EQ(1, pbs_of(result)->n);
    EXPECT_EQ(std::vector<int>({0, 4}), pbs_of(result)->phase_blocks);
    EXPECT_EQ(std::vector<int>({5, 5, 5, 5}), qvars_of(result)->phase_sets);
}

TEST(PhaseblockDataCtor, BoundariesIgnoreUnphasedMiddle) {
    GlobalsGuard guard;
    TempDir dir;

    // a variant with no phase set of its own sits inside the surrounding phase set, so it must not
    // split one block into three: without the backfill running first, every such variant would
    // open a spurious block and collapse the reported phase block sizes
    pipeline_result result = run_pipeline(dir,
            make_qvars({PHASE_ORIG, PHASE_ORIG, PHASE_ORIG}, {7, 0, 7}));
    EXPECT_EQ(1, pbs_of(result)->n);
    EXPECT_EQ(std::vector<int>({0, 3}), pbs_of(result)->phase_blocks);
}

/* write_summary_vcf(): genotype rendering ********************************************************/

/**
 * @brief Builds query variants SPACING bases apart carrying the given per-variant ploidies.
 *
 * A haploid record parses to GT_ALT1_REF on HAP1 alone, exactly as a heterozygous diploid call
 * does, so the ploidy is the only thing distinguishing the two by the time the writer sees them.
 * calc_gts match orig_gts so that every variant classifies as PHASE_ORIG.
 * @param[in] ploidies Ploidy of each variant, in position order
 * @param[in] ctg Contig the variants sit on
 * @return Query variants with orig_gts, calc_gts, phase_sets, and ploidies set
 */
std::shared_ptr<ctgVariants> make_ploidy_qvars(const std::vector<uint8_t> & ploidies,
        const std::string & ctg = CTG) {
    std::vector<var_desc> descs;
    for (size_t i = 0; i < ploidies.size(); i++) {
        var_desc desc;
        desc.pos = int(i) * SPACING;
        desc.rlen = 1;
        desc.ref = "A";
        desc.alt = "C";
        desc.gt = GT_ALT1_REF;
        desc.phase_set = 1;
        desc.ploidy = ploidies[i];
        descs.push_back(desc);
    }
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants(ctg, descs);
    for (size_t i = 0; i < ploidies.size(); i++) qvars->calc_gts[i] = GT_ALT1_REF;
    return qvars;
}

/** @brief Writes the summary VCF for a constructed phaseblockData and returns its contents. */
std::string summary_vcf(const TempDir & dir, phaseblockData & data) {
    std::string vcf_fn = dir.path("summary.vcf");
    std::string log_fn = dir.path("summary.log");
    {
        StderrToFile redirect(log_fn);
        data.write_summary_vcf(vcf_fn);
    }
    return read_text(vcf_fn);
}

/**
 * @brief Returns the query sample's GT field for the record at the given 1-based VCF position.
 * @param[in] vcf Full summary VCF contents
 * @param[in] vcf_pos 1-based POS of the record to read
 * @return The GT field, or the empty string if no record sits at that position
 */
std::string query_gt(const std::string & vcf, int vcf_pos) {
    std::istringstream records(vcf);
    std::string line;
    const std::string prefix = CTG + "\t" + std::to_string(vcf_pos) + "\t";
    while (std::getline(records, line)) {
        if (line.rfind(prefix, 0) != 0) continue;
        size_t sample = line.rfind('\t') + 1;
        return line.substr(sample, line.find(':', sample) - sample);
    }
    return "";
}

TEST(WriteSummaryVcf, HaploidVariantRendersBareAllele) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_ploidy_qvars({1}));
    EXPECT_EQ("1", query_gt(summary_vcf(dir, *result.data), 1));
}

TEST(WriteSummaryVcf, DiploidVariantRendersPhasedPair) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_ploidy_qvars({2}));
    EXPECT_EQ("1|0", query_gt(summary_vcf(dir, *result.data), 1));
}

// A chrX carrying both a PAR diploid call and a non-PAR haploid one renders each in its own shape.
// Both orderings are asserted because per-contig ploidy was inferred from whichever record came
// first, so a single ordering would have agreed with the old behavior half the time.
TEST(WriteSummaryVcf, MixedPloidyHaploidFirst) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_ploidy_qvars({1, 2}));
    std::string vcf = summary_vcf(dir, *result.data);
    EXPECT_EQ("1", query_gt(vcf, 1));
    EXPECT_EQ("1|0", query_gt(vcf, SPACING + 1));
}

TEST(WriteSummaryVcf, MixedPloidyDiploidFirst) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_ploidy_qvars({2, 1}));
    std::string vcf = summary_vcf(dir, *result.data);
    EXPECT_EQ("1|0", query_gt(vcf, 1));
    EXPECT_EQ("1", query_gt(vcf, SPACING + 1));
}

// ploidy= is not a VCF-spec contig attribute, and each record's GT now carries its own.
TEST(WriteSummaryVcf, ContigLineOmitsPloidy) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_ploidy_qvars({2}));
    std::string vcf = summary_vcf(dir, *result.data);
    EXPECT_NE(std::string::npos, vcf.find("##contig=<ID=chr1,length=604>")) << vcf;
    EXPECT_EQ(std::string::npos, vcf.find("ploidy=")) << vcf;
}

} // namespace
