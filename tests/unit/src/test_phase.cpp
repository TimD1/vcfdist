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
 * phase() derives phases[i] from the orig/matched genotype pair rather than reading it, so a test
 * asks for a phasing pattern and gets the genotypes that produce it: 1|0 against 1|0 for
 * PHASE_ORIG, 1|0 against 0|1 for PHASE_SWAP, and 1|1 against 1|1 for PHASE_NONE.
 * @param[in] phases Desired phasing (PHASE_ORIG, PHASE_SWAP, or PHASE_NONE) of each variant
 * @param[in] phase_sets Phase set of each variant, or empty to place them all in phase set 1
 * @param[in] ctg Contig the variants sit on
 * @return Query variants with orig_gts, matched_gts, and phase_sets set
 */
std::shared_ptr<ctgVariants> make_qvars(const std::vector<phase_t> & phases,
        const std::vector<int> & phase_sets = {}, const std::string & ctg = CTG) {
    std::vector<var_desc> descs;
    for (size_t i = 0; i < phases.size(); i++) {
        int phase_set = phase_sets.empty() ? 1 : phase_sets[i];
        gt_t orig_gt = (phases[i] == PHASE_NONE) ? GT_ALT_ALT : GT_ALT_REF;
        descs.push_back({int(i) * SPACING, 1, TYPE_SUB, "A", "C", orig_gt, 60, phase_set});
    }
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants(ctg, descs);
    for (size_t i = 0; i < phases.size(); i++) {
        switch (phases[i]) {
            case PHASE_ORIG: qvars->matched_gts[i] = GT_ALT_REF;  break;
            case PHASE_SWAP: qvars->matched_gts[i] = GT_REF_ALT;  break;
            default:         qvars->matched_gts[i] = GT_ALT_ALT; break;
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
 * @param[in] ref Reference the writer reads anchor bases from, needed only for INS/DEL records
 * @return Constructed phase block data and its captured log
 */
pipeline_result run_pipeline(const TempDir & dir, const std::vector<ctg_input> & inputs,
        std::shared_ptr<fastaData> ref = nullptr) {
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
                make_superclusterData(contigs, lengths, superclusters, ref)));
    }
    result.log = read_text(log_fn);
    return result;
}

/** @brief Runs the phaseblockData pipeline over a single contig. */
pipeline_result run_pipeline(const TempDir & dir, std::shared_ptr<ctgVariants> qvars,
        std::shared_ptr<ctgVariants> tvars = nullptr, int length = CTG_LENGTH,
        std::shared_ptr<fastaData> ref = nullptr) {
    ctg_input input;
    input.qvars = qvars;
    input.tvars = tvars;
    input.length = length;
    return run_pipeline(dir, {input}, ref);
}

/**
 * @brief Builds query variants whose middle variant carries the given genotypes and credits.
 *
 * fix_allele_counts() tie-breaks on the enclosing phase block's phasing, so the middle variant is
 * flanked by two variants that pin the block to block_phase. Every genotype pair passed here
 * classifies as PHASE_NONE, leaving the flanks in sole control of the DP, except a heterozygous
 * call that agrees with itself: that pair is PHASE_ORIG, and is only used in an unswapped block.
 * @param[in] orig_gt Original genotype of the middle variant
 * @param[in] matched_gt Matched genotype of the middle variant
 * @param[in] block_phase Phasing (PHASE_ORIG or PHASE_SWAP) to force on the enclosing block
 * @param[in] hap1_credit Credit of the middle variant on HAP1
 * @param[in] hap2_credit Credit of the middle variant on HAP2
 * @return Query variants whose middle variant is at index 1
 */
std::shared_ptr<ctgVariants> make_ac_qvars(gt_t orig_gt, gt_t matched_gt, phase_t block_phase,
        float hap1_credit = 0, float hap2_credit = 0) {
    std::shared_ptr<ctgVariants> qvars =
            make_qvars({block_phase, PHASE_NONE, block_phase});
    qvars->orig_gts[1] = orig_gt;
    qvars->matched_gts[1] = matched_gt;
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
 * A truth variant's matched_gt is the query genotype recovered for it by alignment, which matched_gts
 * carries in place of the reference call it is initialized to.
 * @param[in] gts Original genotype and recovered query genotype of each variant, in order
 * @return Truth variants with orig_gts and matched_gts set
 */
std::shared_ptr<ctgVariants> make_tvars(const std::vector< std::pair<gt_t, gt_t> > & gts) {
    std::vector<var_desc> descs;
    for (size_t i = 0; i < gts.size(); i++)
        descs.push_back({int(i) * SPACING, 1, TYPE_SUB, "A", "C", gts[i].first, 60, 1});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG, descs);
    for (size_t i = 0; i < gts.size(); i++) tvars->matched_gts[i] = gts[i].second;
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
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 0},
             {100, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 0},
             {200, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 3}});
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
    qvars->matched_gts[0] = GT_REF_REF;
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {CTG}, {CTG_LENGTH}, {make_ctgSuperclusters(qvars, make_ctgVariants(CTG, {}))});
    EXPECT_EXIT(phaseblockData data(sc_data), testing::ExitedWithCode(1),
            "Unknown variant allele count");
}

TEST(FixAlleleCounts, OneToOneTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_REF_ALT, PHASE_ORIG));
    EXPECT_EQ(AC_ERR_1_TO_1, qvars_of(result)->ac_errtype[1]);
    EXPECT_EQ(GT_REF_ALT, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, OneToTwoTallied) {
    GlobalsGuard guard;
    TempDir dir;

    // a heterozygous truth allele called homozygous: one allele is right, the other spurious
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_ALT_ALT, GT_REF_ALT, PHASE_ORIG));
    EXPECT_EQ(AC_ERR_1_TO_2, qvars_of(result)->ac_errtype[1]);
}

TEST(FixAlleleCounts, ZeroToTwoTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_ALT_ALT, GT_REF_REF, PHASE_ORIG));
    EXPECT_EQ(AC_ERR_0_TO_2, qvars_of(result)->ac_errtype[1]);
}

// The truth callset gets the same value the query callset does, read from the other side: a truth
// record's own orig_gt supplies the truth allele count and its recovered matched_gt the query's.

TEST(FixAlleleCounts, TruthHeterozygousFalseNegativeTallied) {
    GlobalsGuard guard;
    TempDir dir;

    // one truth allele that no query allele matched, so its recovered genotype stayed 0|0
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT, GT_REF_REF}}));
    EXPECT_EQ(AC_ERR_1_TO_0, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, TruthHomozygousFalseNegativeTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_ALT_ALT, GT_REF_REF}}));
    EXPECT_EQ(AC_ERR_2_TO_0, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, TruthGenotypeErrorTallied) {
    GlobalsGuard guard;
    TempDir dir;

    // one truth allele the query called on both haplotypes: the truth record reports it too
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT, GT_ALT_ALT}}));
    EXPECT_EQ(AC_ERR_1_TO_2, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, TruthHomozygousMatchTallied) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_ALT_ALT, GT_ALT_ALT}}));
    EXPECT_EQ(AC_ERR_2_TO_2, tvars_of(result)->ac_errtype[0]);
}

TEST(FixAlleleCounts, NoTruthVariantLeftUnknown) {
    GlobalsGuard guard;
    TempDir dir;

    // a pure false negative, a genotype error, and a homozygous match all reach a defined value
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT, GT_REF_REF},
                        {GT_REF_ALT, GT_ALT_ALT},
                        {GT_ALT_ALT, GT_ALT_ALT}}));
    for (int vi = 0; vi < tvars_of(result)->n; vi++)
        EXPECT_NE(AC_UNKNOWN, tvars_of(result)->ac_errtype[vi]) << "truth variant " << vi;
}

TEST(FixAlleleCounts, TruthValuesStayOutOfTheGenotypeErrorSummary) {
    GlobalsGuard guard;
    TempDir dir;

    // the summary tallies the query loop plus the hand-rolled truth false-negative branches, so a
    // truth record now reporting 0/1 -> 1/1 must not add a second count to that row
    pipeline_result result = run_pipeline(dir, nullptr,
            make_tvars({{GT_REF_ALT, GT_ALT_ALT}}));
    EXPECT_TRUE(logged(result.log, "0/1 -> 1/1: 0")) << result.log;
}

TEST(FixAlleleCounts, ForceOneOneKeepsGt) {
    GlobalsGuard guard;
    TempDir dir;

    // a 1|1 call is always evaluated as 1|1, whatever the alignment calculated
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_ALT_ALT, GT_ALT_REF, PHASE_ORIG));
    EXPECT_EQ(GT_ALT_ALT, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, ForceOneOneSwapsHapData) {
    GlobalsGuard guard;
    TempDir dir;

    // a matched_gt of 1|1 carries no record of which haplotype its data came from, so inside a
    // swapped block the per-haplotype lanes are exchanged when the genotype is forced back
    std::shared_ptr<ctgVariants> qvars = make_ac_qvars(GT_ALT_ALT, GT_ALT_ALT, PHASE_SWAP);
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
            make_ac_qvars(GT_REF_ALT, GT_ALT_ALT, PHASE_ORIG, 0.9, 0.1));
    EXPECT_EQ(AC_ERR_2_TO_1, qvars_of(result)->ac_errtype[1]);
    EXPECT_EQ(GT_ALT_REF, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, TwoToOneHap2Better) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_ALT_ALT, PHASE_ORIG, 0.1, 0.9));
    EXPECT_EQ(GT_REF_ALT, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, TwoToOneTieOrig) {
    GlobalsGuard guard;
    TempDir dir;

    // equal credit on both haplotypes falls back to the block's phasing, which keeps the call
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_ALT_ALT, PHASE_ORIG, 0.5, 0.5));
    EXPECT_EQ(GT_REF_ALT, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, TwoToOneTieSwap) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_ALT_ALT, PHASE_SWAP, 0.5, 0.5));
    ASSERT_EQ(PHASE_SWAP, qvars_of(result)->pb_phases[1]);
    EXPECT_EQ(GT_ALT_REF, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneHap1Better) {
    GlobalsGuard guard;
    TempDir dir;

    // a heterozygous call evaluated as reference is restored on the better-credit haplotype
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_REF_REF, PHASE_ORIG, 0.9, 0.1));
    EXPECT_EQ(AC_ERR_0_TO_1, qvars_of(result)->ac_errtype[1]);
    EXPECT_EQ(GT_ALT_REF, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneHap2Better) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_REF_REF, PHASE_ORIG, 0.1, 0.9));
    EXPECT_EQ(GT_REF_ALT, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneTieOrig) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_REF_REF, PHASE_ORIG, 0.5, 0.5));
    EXPECT_EQ(GT_REF_ALT, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, ZeroToOneTieSwap) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir,
            make_ac_qvars(GT_REF_ALT, GT_REF_REF, PHASE_SWAP, 0.5, 0.5));
    ASSERT_EQ(PHASE_SWAP, qvars_of(result)->pb_phases[1]);
    EXPECT_EQ(GT_ALT_REF, qvars_of(result)->matched_gts[1]);
}

TEST(FixAlleleCounts, TruthFn2To0) {
    GlobalsGuard guard;
    TempDir dir;

    // a homozygous truth variant missed on both haplotypes loses two alleles, and only the truth
    // VCF can show it: there is no query variant to tally it against
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG,
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT_ALT, 60, 1}});
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
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 1}});
    set_hap_data(tvars, HAP1, 0, ERRTYPE_FN, 0, 0, 1, 1, 0);
    pipeline_result result = run_pipeline(dir, nullptr, tvars);
    EXPECT_TRUE(logged(result.log, "0/1 -> 0/0: 1 "));
}

TEST(FixAlleleCounts, TruthFn1To0Hap2) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants(CTG,
            {{0, 1, TYPE_SUB, "A", "C", GT_REF_ALT, 60, 1}});
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
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 1},
             {100, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 1},
             {200, 1, TYPE_SUB, "A", "C", GT_ALT_ALT, 60, 1}});
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
            {{0, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 1},
             {100, 1, TYPE_SUB, "A", "C", GT_ALT_REF, 60, 1}});
    second.qvars->matched_gts[0] = GT_ALT_REF;
    second.qvars->matched_gts[1] = GT_ALT_REF;
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
 * A haploid record parses to GT_ALT_REF on HAP1 alone, exactly as a heterozygous diploid call
 * does, so the ploidy is the only thing distinguishing the two by the time the writer sees them.
 * matched_gts match orig_gts so that every variant classifies as PHASE_ORIG.
 * @param[in] ploidies Ploidy of each variant, in position order
 * @param[in] ctg Contig the variants sit on
 * @return Query variants with orig_gts, matched_gts, phase_sets, and ploidies set
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
        desc.gt = GT_ALT_REF;
        desc.phase_set = 1;
        desc.ploidy = ploidies[i];
        descs.push_back(desc);
    }
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants(ctg, descs);
    for (size_t i = 0; i < ploidies.size(); i++) qvars->matched_gts[i] = GT_ALT_REF;
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

/* write_summary_vcf(): one record per variant *****************************************************/

const int QUERY_COL = 10; ///< 0-based column of the QUERY sample within a summary VCF record

/** @brief Field offsets within a summary VCF sample column, in FORMAT order. */
enum sample_field { FMT_GT, FMT_BD, FMT_BC, FMT_RD, FMT_QD, FMT_BK, FMT_QQ, FMT_SC, FMT_SG };

/** @brief Returns every summary VCF record at the given 1-based position, in file order. */
std::vector<std::string> records_at(const std::string & vcf, int vcf_pos) {
    std::vector<std::string> found;
    std::istringstream lines(vcf);
    std::string line;
    const std::string prefix = CTG + "\t" + std::to_string(vcf_pos) + "\t";
    while (std::getline(lines, line)) {
        if (line.rfind(prefix, 0) == 0) found.push_back(line);
    }
    return found;
}

/** @brief Splits a whole record on a delimiter, yielding "" for a column past the end. */
std::vector<std::string> split(const std::string & text, char delim) {
    std::vector<std::string> parts;
    std::istringstream fields(text);
    std::string part;
    while (std::getline(fields, part, delim)) parts.push_back(part);
    return parts;
}

/**
 * @brief Returns the QUERY sample's FORMAT fields for the sole record at a 1-based position.
 *
 * The record count is asserted here rather than returned, so that a test reading a field also
 * pins that exactly one record was written for the variant.
 * @param[in] vcf Full summary VCF contents
 * @param[in] vcf_pos 1-based POS of the record to read
 * @return The QUERY sample column split on ':', indexable by sample_field
 */
std::vector<std::string> sole_query_sample(const std::string & vcf, int vcf_pos = 1) {
    std::vector<std::string> recs = records_at(vcf, vcf_pos);
    EXPECT_EQ(size_t(1), recs.size()) << vcf;
    if (recs.size() != 1) return std::vector<std::string>(15, "");
    return split(split(recs[0], '\t').at(QUERY_COL), ':');
}

/**
 * @brief Builds one query variant of the given type and genotypes, SPACING bases into the contig.
 *
 * An INS/DEL is left-anchored on the preceding reference base, so it cannot sit at position 0 the
 * way the substitution builders above place their first variant.
 * @param[in] type Variant type (TYPE_SUB, TYPE_INS, or TYPE_DEL)
 * @param[in] orig_gt Original genotype, reported in the record's GT column
 * @param[in] matched_gt Calculated genotype, which indexes the per-haplotype evaluation lanes
 * @param[in] ploidy Variant ploidy (0 = unknown, treated as diploid)
 * @return Query variants holding the single described variant
 */
std::shared_ptr<ctgVariants> make_shape_qvars(edittype_t type, gt_t orig_gt, gt_t matched_gt,
        uint8_t ploidy = 2) {
    var_desc desc;
    desc.pos = SPACING;
    desc.rlen = type == TYPE_INS ? 0 : 1;
    desc.type = type;
    desc.ref = type == TYPE_INS ? "" : "A";
    desc.alt = type == TYPE_DEL ? "" : "C";
    desc.gt = orig_gt;
    desc.phase_set = 1;
    desc.ploidy = ploidy;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants(CTG, {desc});
    qvars->matched_gts[0] = matched_gt;
    return qvars;
}

/** @brief Runs the pipeline over one contig of all-'A' reference and returns the summary VCF. */
std::string shape_vcf(const TempDir & dir, std::shared_ptr<ctgVariants> qvars) {
    pipeline_result result = run_pipeline(dir, qvars, nullptr, CTG_LENGTH,
            make_fasta(CTG, std::string(CTG_LENGTH, 'A')));
    return summary_vcf(dir, *result.data);
}

// A homozygous variant reaches the writer as a single entry carrying GT_ALT_ALT, so writing one
// record per haplotype re-split it into two. Both alleles are alternate, so both carry data.
TEST(WriteSummaryVcf, HomSnpIsOneRecordWithTwoValues) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars =
            make_shape_qvars(TYPE_SUB, GT_ALT_ALT, GT_ALT_ALT);
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 4, 60, 3, 0, 1.0);
    set_hap_data(qvars, HAP2, 0, ERRTYPE_TP, 4, 60, 3, 0, 1.0);

    std::vector<std::string> sample = sole_query_sample(shape_vcf(dir, qvars), SPACING + 1);
    EXPECT_EQ("1|1", sample.at(FMT_GT));
    EXPECT_EQ("TP,TP", sample.at(FMT_BD));
    EXPECT_EQ("1.000000,1.000000", sample.at(FMT_BC));
    EXPECT_EQ("3,3", sample.at(FMT_RD));
    EXPECT_EQ("0,0", sample.at(FMT_QD));
    EXPECT_EQ("gm,gm", sample.at(FMT_BK));
    EXPECT_EQ("4,4", sample.at(FMT_SG));
}

// An indel takes the same shape as a substitution; it is written through print_var_info's
// left-anchoring branch, so it is covered separately.
TEST(WriteSummaryVcf, HomIndelIsOneRecordWithTwoValues) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars =
            make_shape_qvars(TYPE_DEL, GT_ALT_ALT, GT_ALT_ALT);
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 0, 60, 1, 0, 1.0);
    set_hap_data(qvars, HAP2, 0, ERRTYPE_TP, 0, 60, 1, 0, 1.0);

    // a deletion is positioned on its anchor base, one before its first deleted base
    std::vector<std::string> recs = records_at(shape_vcf(dir, qvars), SPACING);
    ASSERT_EQ(size_t(1), recs.size());
    std::vector<std::string> cols = split(recs[0], '\t');
    EXPECT_EQ("AA", cols.at(3));
    EXPECT_EQ("A", cols.at(4));
    std::vector<std::string> sample = split(cols.at(QUERY_COL), ':');
    EXPECT_EQ("1|1", sample.at(FMT_GT));
    EXPECT_EQ("TP,TP", sample.at(FMT_BD));
    EXPECT_EQ("1.000000,1.000000", sample.at(FMT_BC));
}

// A heterozygous call still carries one value per GT allele, but its reference allele was never
// evaluated, so reporting the untouched haplotype's lane there would invent a second decision.
TEST(WriteSummaryVcf, HetRecordDotsTheReferenceAllele) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars =
            make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF);
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 2, 60, 5, 0, 1.0);

    std::vector<std::string> sample = sole_query_sample(shape_vcf(dir, qvars), SPACING + 1);
    EXPECT_EQ("1|0", sample.at(FMT_GT));
    EXPECT_EQ("TP,.", sample.at(FMT_BD));
    EXPECT_EQ("1.000000,.", sample.at(FMT_BC));
    EXPECT_EQ("5,.", sample.at(FMT_RD));
    EXPECT_EQ("2,.", sample.at(FMT_SG));
}

// A haploid record has one GT allele, so it carries one value rather than a two-value list.
TEST(WriteSummaryVcf, HaploidRecordCarriesOneValue) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars =
            make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF, 1);
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 2, 60, 5, 0, 1.0);

    std::vector<std::string> sample = sole_query_sample(shape_vcf(dir, qvars), SPACING + 1);
    EXPECT_EQ("1", sample.at(FMT_GT));
    EXPECT_EQ("TP", sample.at(FMT_BD));
    EXPECT_EQ("1.000000", sample.at(FMT_BC));
    EXPECT_EQ("5", sample.at(FMT_RD));
    EXPECT_EQ("2", sample.at(FMT_SG));
}

// The values are ordered by GT allele, and GT reports orig_gt while the evaluation lanes are
// indexed by matched_gt's haplotypes. On a swapped record the two disagree, so reading the lanes in
// haplotype order would report the untouched haplotype's zero credit against the alternate allele.
TEST(WriteSummaryVcf, PerAlleleValuesFollowTheGenotypeSwap) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars =
            make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_REF_ALT);
    set_hap_data(qvars, HAP1, 0, ERRTYPE_FP, 7, 60, 0, 0, 0.0);
    set_hap_data(qvars, HAP2, 0, ERRTYPE_TP, 2, 60, 5, 0, 1.0);

    std::vector<std::string> sample = sole_query_sample(shape_vcf(dir, qvars), SPACING + 1);
    EXPECT_EQ("1|0", sample.at(FMT_GT));
    EXPECT_EQ("TP,.", sample.at(FMT_BD));
    EXPECT_EQ("1.000000,.", sample.at(FMT_BC));
    EXPECT_EQ("2,.", sample.at(FMT_SG));
}

// Both alleles of a homozygous call carry data, so their order is directly observable.
TEST(WriteSummaryVcf, PerAlleleValuesAreInGenotypeAlleleOrder) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars =
            make_shape_qvars(TYPE_SUB, GT_ALT_ALT, GT_ALT_ALT);
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 1, 60, 3, 0, 1.0);
    set_hap_data(qvars, HAP2, 0, ERRTYPE_TP, 6, 60, 7, 2, 0.75);

    std::vector<std::string> sample = sole_query_sample(shape_vcf(dir, qvars), SPACING + 1);
    EXPECT_EQ("3,7", sample.at(FMT_RD));
    EXPECT_EQ("0,2", sample.at(FMT_QD));
    EXPECT_EQ("1,6", sample.at(FMT_SG));
    EXPECT_EQ("1.000000,0.750000", sample.at(FMT_BC));
}

/**
 * @brief Builds the two co-located entries a het-alt (1|2) source record is parsed into.
 *
 * Both entries share one record ordinal but carry different ALT ordinals, which is what an
 * ALT-indexed source value has to be subset against.
 * @return Query variants holding the entries for ALT ordinals 1 and 2, in that order
 */
std::shared_ptr<ctgVariants> make_het_alt_qvars() {
    std::vector<var_desc> descs;
    for (const std::string & alt : {"C", "G"}) {
        var_desc desc;
        desc.pos = SPACING;
        desc.rlen = 1;
        desc.ref = "A";
        desc.alt = alt;
        desc.gt = alt == "C" ? GT_ALT_REF : GT_REF_ALT;
        desc.phase_set = 1;
        desc.alt_idx = alt == "C" ? 1 : 2;
        desc.ploidy = 2;
        descs.push_back(desc);
    }
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants(CTG, descs);
    qvars->matched_gts[0] = GT_ALT_REF;
    qvars->matched_gts[1] = GT_REF_ALT;
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 0, 60, 1, 0, 1.0);
    set_hap_data(qvars, HAP2, 1, ERRTYPE_TP, 0, 60, 1, 0, 1.0);
    return qvars;
}

// A het-alt (1|2) source record is parsed into two entries with different ALTs, and the
// cross-haplotype merge only collapses entries whose position and alleles match exactly, so
// nothing rejoins them. Pinned explicitly: this asymmetry against the homozygous case is
// deliberate, since the two alleles need not even share a position once normalized.
TEST(WriteSummaryVcf, HetAltStaysTwoColocatedRecords) {
    GlobalsGuard guard;
    TempDir dir;
    std::vector<std::string> recs = records_at(shape_vcf(dir, make_het_alt_qvars()), SPACING + 1);
    ASSERT_EQ(size_t(2), recs.size());
    EXPECT_EQ("C", split(recs[0], '\t').at(4));
    EXPECT_EQ("G", split(recs[1], '\t').at(4));
    EXPECT_EQ("1|0", split(split(recs[0], '\t').at(QUERY_COL), ':').at(FMT_GT));
    EXPECT_EQ("0|1", split(split(recs[1], '\t').at(QUERY_COL), ':').at(FMT_GT));
}

// Number=P would declare the one-value-per-GT-allele cardinality these fields carry, but it is a
// VCF 4.4 addition htslib only supports from 1.23, so Number=. is declared and the count and
// order are stated in the description instead.
TEST(WriteSummaryVcf, PerAlleleFieldsAreDeclaredUnbounded) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = shape_vcf(dir, make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF));
    for (const std::string & id : {"BD", "BC", "RD", "QD", "BK", "SG"}) {
        EXPECT_NE(std::string::npos, vcf.find("##FORMAT=<ID=" + id + ",Number=.,")) << id;
        EXPECT_NE(std::string::npos,
                vcf.find("One value per allele of this sample's GT, in GT allele order")) << id;
    }
    for (const std::string & id : {"GT", "QQ", "SC", "PS", "PB", "BS", "VP", "FE", "GE"}) {
        EXPECT_NE(std::string::npos, vcf.find("##FORMAT=<ID=" + id + ",Number=1,")) << id;
    }
}

/* write_summary_vcf(): preserved source fields ***************************************************/

const int ID_COL = 2;     ///< 0-based column of ID within a summary VCF record
const int ALT_COL = 4;    ///< 0-based column of ALT within a summary VCF record
const int QUAL_COL = 5;   ///< 0-based column of QUAL within a summary VCF record
const int FILTER_COL = 6; ///< 0-based column of FILTER within a summary VCF record
const int INFO_COL = 7;   ///< 0-based column of INFO within a summary VCF record
const int FORMAT_COL = 8; ///< 0-based column of FORMAT within a summary VCF record
const int TRUTH_COL = 9;  ///< 0-based column of the TRUTH sample within a summary VCF record

/**
 * @brief Attaches a one-record source store to a variant container, at record ordinal 0.
 * @param[in,out] vars Container whose sole variant is pointed at the retained record
 * @param[in] id ID column of the retained record
 * @param[in] info INFO column of the retained record
 * @param[in] keys Preserved FORMAT keys, each prefixed with ':'
 * @param[in] vals Preserved FORMAT values, each prefixed with ':'
 * @return The store, so a caller can also hand it to phaseblockData for its header lines
 */
std::shared_ptr<srcRecords> attach_src(std::shared_ptr<ctgVariants> vars, const std::string & id,
        const std::string & info, const std::string & keys, const std::string & vals) {
    vars->src_recs = std::shared_ptr<srcRecords>(new srcRecords());
    vars->src_recs->add(0, id, "37", "LowConf", info, keys, vals);
    for (int i = 0; i < vars->n; i++) vars->rec_idxs[i] = 0;
    return vars->src_recs;
}

/** @brief Builds one truth variant matching make_shape_qvars()'s substitution, allele for allele. */
std::shared_ptr<ctgVariants> make_shape_tvars(gt_t orig_gt) {
    var_desc desc;
    desc.pos = SPACING;
    desc.rlen = 1;
    desc.ref = "A";
    desc.alt = "C";
    desc.gt = orig_gt;
    desc.phase_set = 1;
    desc.ploidy = 2;
    return make_ctgVariants(CTG, {desc});
}

// Before any field was preserved these four columns were hardcoded, so a variant whose source
// record was never retained must still render exactly what it used to.
TEST(WriteSummaryVcf, PlaceholdersWhenNoSourceRecordWasRetained) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = shape_vcf(dir, make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF));
    std::vector<std::string> cols = split(records_at(vcf, SPACING + 1).at(0), '\t');
    EXPECT_EQ(".", cols.at(ID_COL));
    EXPECT_EQ(".", cols.at(QUAL_COL));
    EXPECT_EQ("PASS", cols.at(FILTER_COL));
    EXPECT_EQ(".", cols.at(INFO_COL));
    EXPECT_EQ("GT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE", cols.at(FORMAT_COL));
}

TEST(WriteSummaryVcf, SourceColumnsReplacePlaceholders) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars = make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF);
    attach_src(qvars, "rs1", "DP=30;SOMATIC", ":SDP:SAC", ":29:12,17");

    std::vector<std::string> cols = split(records_at(shape_vcf(dir, qvars),
            SPACING + 1).at(0), '\t');
    EXPECT_EQ("rs1", cols.at(ID_COL));
    EXPECT_EQ("37", cols.at(QUAL_COL));
    EXPECT_EQ("LowConf", cols.at(FILTER_COL));
    EXPECT_EQ("DP=30;SOMATIC", cols.at(INFO_COL));
    EXPECT_EQ("GT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE:SDP:SAC", cols.at(FORMAT_COL));

    // the appended values close out the query sample, in the order the keys were appended
    std::vector<std::string> sample = split(cols.at(QUERY_COL), ':');
    EXPECT_EQ("29", sample.at(sample.size() - 2));
    EXPECT_EQ("12,17", sample.at(sample.size() - 1));
}

// The two entries of a het-alt record share one retained copy of the source columns, so an
// ALT-indexed value only becomes correct once it is subset: each output line must receive the
// elements of its own ALT rather than the whole list the source record wrote.
TEST(WriteSummaryVcf, HetAltEntriesReceiveDifferentAltIndexedValues) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars = make_het_alt_qvars();
    std::shared_ptr<srcRecords> src = attach_src(qvars, "rs250", "DP=32;AF=0.4,0.6",
            ":AD:PL", ":0,15,16:255,60,0,44,11,7");
    src->info_lens["AF"] = BCF_VL_A;
    src->fmt_lens["AD"] = BCF_VL_R;
    src->fmt_lens["PL"] = BCF_VL_G;

    std::vector<std::string> recs = records_at(shape_vcf(dir, qvars), SPACING + 1);
    ASSERT_EQ(size_t(2), recs.size());
    std::vector<std::string> first = split(recs[0], '\t');
    std::vector<std::string> second = split(recs[1], '\t');
    EXPECT_EQ("DP=32;AF=0.4", first.at(INFO_COL));
    EXPECT_EQ("DP=32;AF=0.6", second.at(INFO_COL));

    // AD and PL close out the sample column, in the order their keys were appended
    std::vector<std::string> first_sample = split(first.at(QUERY_COL), ':');
    std::vector<std::string> second_sample = split(second.at(QUERY_COL), ':');
    EXPECT_EQ("0,15", first_sample.at(first_sample.size() - 2));
    EXPECT_EQ("255,60,0", first_sample.back());
    EXPECT_EQ("0,16", second_sample.at(second_sample.size() - 2));
    EXPECT_EQ("255,44,7", second_sample.back());
}

// One FORMAT key list serves both samples, so the callset that did not supply the appended keys
// has no values for them and reports each as missing.
TEST(WriteSummaryVcf, TruthSamplePadsQueryOwnedFormatKeys) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars = make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF);
    attach_src(qvars, "rs1", "DP=30", ":SDP:SAC", ":29:12,17");

    pipeline_result result = run_pipeline(dir, qvars, make_shape_tvars(GT_ALT_REF), CTG_LENGTH,
            make_fasta(CTG, std::string(CTG_LENGTH, 'A')));
    std::vector<std::string> cols = split(records_at(summary_vcf(dir, *result.data),
            SPACING + 1).at(0), '\t');
    std::vector<std::string> truth = split(cols.at(TRUTH_COL), ':');
    EXPECT_EQ(".", truth.at(truth.size() - 2));
    EXPECT_EQ(".", truth.at(truth.size() - 1));
}

// A record the query never calls on is owned by the truth, so its columns and appended keys are
// the truth's and the query sample is the one that pads.
TEST(WriteSummaryVcf, TruthOnlyRecordCarriesTruthColumns) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> tvars = make_shape_tvars(GT_ALT_REF);
    attach_src(tvars, "tv1", "TRUTHSET=hc", ":SDP", ":98");

    pipeline_result result = run_pipeline(dir, nullptr, tvars, CTG_LENGTH,
            make_fasta(CTG, std::string(CTG_LENGTH, 'A')));
    std::vector<std::string> cols = split(records_at(summary_vcf(dir, *result.data),
            SPACING + 1).at(0), '\t');
    EXPECT_EQ("tv1", cols.at(ID_COL));
    EXPECT_EQ("TRUTHSET=hc", cols.at(INFO_COL));
    EXPECT_EQ("GT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE:SDP", cols.at(FORMAT_COL));
    EXPECT_EQ("98", split(cols.at(TRUTH_COL), ':').back());
    EXPECT_EQ(".", split(cols.at(QUERY_COL), ':').back());
}

// Both callsets contribute declarations, but a header may declare an ID only once; the query's
// wins, matching which callset supplies the columns of a record both call.
TEST(WriteSummaryVcf, HeaderDeclaresEachRetainedFieldOnce) {
    GlobalsGuard guard;
    TempDir dir;
    pipeline_result result = run_pipeline(dir, make_shape_qvars(TYPE_SUB, GT_ALT_REF,
            GT_ALT_REF), nullptr, CTG_LENGTH, make_fasta(CTG, std::string(CTG_LENGTH, 'A')));

    std::shared_ptr<srcRecords> qsrc(new srcRecords());
    qsrc->hdr_keys = {"INFO/DP", "FILTER/LowConf"};
    qsrc->hdr_lines = {"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Query depth\">",
                       "##FILTER=<ID=LowConf,Description=\"Low confidence call\">"};
    std::shared_ptr<srcRecords> tsrc(new srcRecords());
    tsrc->hdr_keys = {"INFO/DP", "INFO/TRUTHSET"};
    tsrc->hdr_lines = {"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Truth depth\">",
                       "##INFO=<ID=TRUTHSET,Number=1,Type=String,Description=\"Truth set\">"};
    result.data->callset_src_recs = {qsrc, tsrc};

    std::string vcf = summary_vcf(dir, *result.data);
    EXPECT_NE(std::string::npos, vcf.find("Query depth")) << vcf;
    EXPECT_EQ(std::string::npos, vcf.find("Truth depth")) << vcf;
    EXPECT_NE(std::string::npos, vcf.find("##FILTER=<ID=LowConf,")) << vcf;
    EXPECT_NE(std::string::npos, vcf.find("##INFO=<ID=TRUTHSET,")) << vcf;

    // the writer declares PASS itself, so an input declaration of it would be a duplicate
    size_t first_pass = vcf.find("##FILTER=<ID=PASS,");
    ASSERT_NE(std::string::npos, first_pass);
    EXPECT_EQ(std::string::npos, vcf.find("##FILTER=<ID=PASS,", first_pass + 1));
}

/* write_summary_vcf(): retained but unevaluated records *******************************************/

/**
 * @brief Builds a sideline container holding one A>C record per (position, reason) pair.
 * @param[in] retained Position (0-based) and SIDELINE_* reason of each retained record, in order
 * @return Retained records of one contig, each keyed by its own record ordinal
 */
std::shared_ptr<ctgSideline> make_sideline(
        const std::vector< std::pair<int, uint8_t> > & retained) {
    std::shared_ptr<ctgSideline> side(new ctgSideline(CTG));
    for (size_t i = 0; i < retained.size(); i++)
        side->add(int(i), SIDELINE_ALL_HAPS, retained[i].first, "A", "C", "1|0",
                retained[i].second);
    return side;
}

/**
 * @brief Runs the pipeline over one query variant plus the given retained records of each callset.
 * @param[in] dir Temporary directory receiving the pipeline's output
 * @param[in] qvars Query variants to evaluate, or nullptr for none
 * @param[in] qside Retained query records, or nullptr for none
 * @param[in] tside Retained truth records, or nullptr for none
 * @return The summary VCF the writer produced
 */
std::string sideline_vcf(const TempDir & dir, std::shared_ptr<ctgVariants> qvars,
        std::shared_ptr<ctgSideline> qside, std::shared_ptr<ctgSideline> tside = nullptr) {
    pipeline_result result = run_pipeline(dir, qvars, nullptr, CTG_LENGTH,
            make_fasta(CTG, std::string(CTG_LENGTH, 'A')));
    std::unordered_map< std::string, std::shared_ptr<ctgSideline> > qmap, tmap;
    if (qside != nullptr) qmap[CTG] = qside;
    if (tside != nullptr) tmap[CTG] = tside;
    result.data->callset_sidelined = {qmap, tmap};
    return summary_vcf(dir, *result.data);
}

// A record excluded from evaluation is still written, reporting BD=N and nothing else: every other
// fixed field is the result of an evaluation that never ran.
TEST(WriteSummaryVcf, RetainedRecordReportsNotAssessed) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = sideline_vcf(dir, nullptr,
            make_sideline({{SPACING, SIDELINE_LOW_QUAL}}));

    std::vector<std::string> cols = split(records_at(vcf, SPACING + 1).at(0), '\t');
    EXPECT_EQ("VCFDIST_LOW_QUAL", cols.at(FILTER_COL));
    std::vector<std::string> query = split(cols.at(QUERY_COL), ':');
    EXPECT_EQ("1|0", query.at(FMT_GT));
    EXPECT_EQ("N", query.at(FMT_BD));
    for (size_t i = FMT_BC; i < query.size(); i++) EXPECT_EQ(".", query.at(i)) << i;

    // it was excluded before any comparison ran, so it can never have been matched
    EXPECT_EQ(std::string(".:.:.:.:.:.:.:.:.:.:.:.:.:.:."), cols.at(TRUTH_COL));
}

// The tag names why the record was not evaluated, and is prefixed so that no input FILTER ID can
// collide with it.
TEST(WriteSummaryVcf, RetainedRecordFilterNamesTheReason) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = sideline_vcf(dir, nullptr,
            make_sideline({{SPACING, SIDELINE_FAILED_FILTER}}));
    EXPECT_EQ("VCFDIST_FAILED_FILTER",
            split(records_at(vcf, SPACING + 1).at(0), '\t').at(FILTER_COL));
}

// FILTER lists the filters a record failed, so the tag joins the record's own list rather than
// replacing it; only a lone PASS, which a retained record contradicts, is replaced.
TEST(WriteSummaryVcf, RetainedRecordKeepsItsInputFilters) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgSideline> side = make_sideline({{SPACING, SIDELINE_LOW_QUAL}});
    side->src_recs = std::shared_ptr<srcRecords>(new srcRecords());
    side->src_recs->add(0, "rs1", "5", "LowConf", "DP=30", ":SDP", ":29");

    std::vector<std::string> cols = split(records_at(sideline_vcf(dir, nullptr, side),
            SPACING + 1).at(0), '\t');
    EXPECT_EQ("LowConf;VCFDIST_LOW_QUAL", cols.at(FILTER_COL));
    EXPECT_EQ("rs1", cols.at(ID_COL));
    EXPECT_EQ("5", cols.at(QUAL_COL));
    EXPECT_EQ("DP=30", cols.at(INFO_COL));
    EXPECT_EQ("GT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE:SDP", cols.at(FORMAT_COL));
    EXPECT_EQ("29", split(cols.at(QUERY_COL), ':').back());
    EXPECT_EQ(".", split(cols.at(TRUTH_COL), ':').back());
}

// Retained records are merged into the evaluated walk by position, not appended to it, so the file
// stays sorted whether a retained record precedes or follows the calls around it.
TEST(WriteSummaryVcf, RetainedRecordsInterleaveByPosition) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = sideline_vcf(dir, make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF),
            make_sideline({{0, SIDELINE_LOW_QUAL}, {2*SPACING, SIDELINE_FAILED_FILTER}}));

    std::vector<int> written;
    std::istringstream lines(vcf);
    std::string line;
    while (std::getline(lines, line)) {
        if (line.rfind(CTG + "\t", 0) != 0) continue;
        written.push_back(std::stoi(split(line, '\t').at(1)));
    }
    EXPECT_EQ(std::vector<int>({1, SPACING + 1, 2*SPACING + 1}), written);
}

// Each callset's retained records are its own, so the sample carrying the call is the one that
// called it and the other stays empty.
TEST(WriteSummaryVcf, TruthRetainedRecordCarriesTheTruthCall) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = sideline_vcf(dir, nullptr, nullptr,
            make_sideline({{SPACING, SIDELINE_LOW_QUAL}}));

    std::vector<std::string> cols = split(records_at(vcf, SPACING + 1).at(0), '\t');
    EXPECT_EQ("1|0", split(cols.at(TRUTH_COL), ':').at(FMT_GT));
    EXPECT_EQ("N", split(cols.at(TRUTH_COL), ':').at(FMT_BD));
    EXPECT_EQ(std::string(".:.:.:.:.:.:.:.:.:.:.:.:.:.:."), cols.at(QUERY_COL));
}

// Every tag a retained record can carry is declared, whether or not this run wrote one, and BD's
// description covers the value they report.
TEST(WriteSummaryVcf, HeaderDeclaresEveryRetentionReason) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = sideline_vcf(dir, make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF),
            nullptr);
    for (const std::string & tag : sideline_strs)
        EXPECT_NE(std::string::npos, vcf.find("##FILTER=<ID=" + tag + ",Description=\"")) << tag;
    EXPECT_NE(std::string::npos, vcf.find("N for a call that was not assessed")) << vcf;
}

// Each reason writes a tag of its own, so a user reading the output can tell which check excluded
// the record; two reasons sharing a tag would make the three BED conditions indistinguishable.
TEST(WriteSummaryVcf, EachRetentionReasonWritesItsOwnTag) {
    GlobalsGuard guard;
    for (uint8_t reason = 0; reason < SIDELINES; reason++) {
        TempDir dir;
        std::vector<std::string> recs = records_at(
                sideline_vcf(dir, nullptr, make_sideline({{SPACING, reason}})), SPACING + 1);
        ASSERT_EQ(size_t(1), recs.size()) << int(reason);
        EXPECT_EQ(sideline_strs[reason], split(recs[0], '\t').at(FILTER_COL)) << int(reason);
    }
}

/**
 * @brief Runs the pipeline over one query variant plus retained records on a second contig.
 * @param[in] dir Temporary directory receiving the pipeline's output
 * @param[in] off Retained query records of a contig the evaluated walk does not cover
 * @return The summary VCF the writer produced
 */
std::string off_contig_vcf(const TempDir & dir, std::shared_ptr<ctgSideline> off) {
    pipeline_result result = run_pipeline(dir,
            make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF), nullptr, CTG_LENGTH,
            make_fasta(CTG, std::string(CTG_LENGTH, 'A')));
    std::unordered_map< std::string, std::shared_ptr<ctgSideline> > qmap;
    qmap[off->ctg] = off;
    result.data->callset_sidelined = {qmap,
            std::unordered_map< std::string, std::shared_ptr<ctgSideline> >()};
    return summary_vcf(dir, *result.data);
}

// A contig absent from the BED file is dropped before evaluation, so the evaluated walk never
// reaches it, yet its records are retained and still belong in the output. They are written after
// the evaluated contigs, and the contig is declared in the same order: a record on an undeclared
// contig is not a valid VCF, and a header ordering the contigs otherwise would unsort the file.
TEST(WriteSummaryVcf, RetainedRecordsOnAnUncoveredContigAreStillWritten) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgSideline> off(new ctgSideline("chr9", 3, 250));
    off->add(0, SIDELINE_ALL_HAPS, 40, "T", "C", "1|1", SIDELINE_BED_OFF_CTG);
    std::string vcf = off_contig_vcf(dir, off);

    EXPECT_NE(std::string::npos, vcf.find("##contig=<ID=chr9,length=250>")) << vcf;
    EXPECT_LT(vcf.find("##contig=<ID=" + CTG + ","), vcf.find("##contig=<ID=chr9,")) << vcf;

    std::string record;
    std::istringstream lines(vcf);
    std::string line;
    while (std::getline(lines, line)) if (line.rfind("chr9\t", 0) == 0) record = line;
    ASSERT_FALSE(record.empty()) << vcf;
    std::vector<std::string> cols = split(record, '\t');
    EXPECT_EQ("41", cols.at(1));
    EXPECT_EQ("VCFDIST_BED_OFF_CTG", cols.at(FILTER_COL));
    EXPECT_EQ("N", split(cols.at(QUERY_COL), ':').at(FMT_BD));
}

// Only a contig that retained something is declared, since a container exists for every contig of
// every input and declaring them all would add contigs the run never had a record on.
TEST(WriteSummaryVcf, AnUncoveredContigRetainingNothingIsNotDeclared) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf = off_contig_vcf(dir, std::shared_ptr<ctgSideline>(
            new ctgSideline("chr9", 3, 250)));
    EXPECT_EQ(std::string::npos, vcf.find("##contig=<ID=chr9,")) << vcf;
}

// A per-allele reason can exclude one haplotype of a record while the other stays evaluated, so
// both land at the record's position. The evaluated record is written first: a retained one
// belongs to no supercluster, so it has no place within the evaluated ordering of a position. It
// carries the whole ALT list and genotype of the record it came from, neither having been split.
TEST(WriteSummaryVcf, PartiallyRetainedRecordIsWrittenBesideItsEvaluatedAllele) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgSideline> side(new ctgSideline(CTG));
    side->add(0, int(HAP2), SPACING, "A", "C,CCCCCCCC", "1|2", SIDELINE_TOO_LARGE);
    std::vector<std::string> recs = records_at(
            sideline_vcf(dir, make_shape_qvars(TYPE_SUB, GT_ALT_REF, GT_ALT_REF), side),
            SPACING + 1);

    ASSERT_EQ(size_t(2), recs.size());
    EXPECT_EQ("PASS", split(recs[0], '\t').at(FILTER_COL));
    EXPECT_EQ("C", split(recs[0], '\t').at(ALT_COL));
    EXPECT_EQ("VCFDIST_TOO_LARGE", split(recs[1], '\t').at(FILTER_COL));
    EXPECT_EQ("C,CCCCCCCC", split(recs[1], '\t').at(ALT_COL));
    EXPECT_EQ("1|2", split(split(recs[1], '\t').at(QUERY_COL), ':').at(FMT_GT));
    EXPECT_EQ("N", split(split(recs[1], '\t').at(QUERY_COL), ':').at(FMT_BD));
}

} // namespace
