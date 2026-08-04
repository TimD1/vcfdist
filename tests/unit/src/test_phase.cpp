/**
 * @file test_phase.cpp
 * @brief Unit tests for phase.cpp: phase block statistics and switch/flip classification.
 */
#include "gtest/gtest.h"

#include "../../../src/globals.h"
#include "../../../src/phase.h"
#include "test_helpers.h"

namespace {

/* correct_block_sizes ****************************************************************************/

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

} // namespace
