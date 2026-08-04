/**
 * @file test_dist.cpp
 * @brief Unit tests for dist.cpp: alignment primitives, NG50, graph cell indices, and the
 *        graph-based precision/recall evaluation path.
 *
 * Every expected edit distance, alignment score, credit and edit-distance figure below is derived
 * by hand from the penalties the algorithm charges (substitution x, gap of length L costs
 * o + e*L for wf_swg_align; unit cost per edit plus a (1 - credit_threshold) * len bypass toll
 * for calc_prec_recall_aln), never by recording what the implementation happens to print. The
 * derivation is spelled out in a comment wherever the arithmetic is not immediate.
 */
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "../../../src/dist.h"
#include "../../../src/globals.h"
#include "test_helpers.h"

namespace {

/** @brief Sets the three Smith-Waterman penalties so no test depends on today's defaults. */
static void set_penalties(int sub, int open, int extend) {
    g.sub = sub;
    g.open = open;
    g.extend = extend;
}

/* contains ***************************************************************************************/

// contains is instantiated below over element types dist.cpp never uses, so these tests only link
// while its definition lives in dist.h rather than dist.cpp.

TEST(Contains, SetInt) {
    std::unordered_set<int> hash_set = {-7, 0, 42};
    EXPECT_TRUE(contains(hash_set, -7));
    EXPECT_TRUE(contains(hash_set, 0));
    EXPECT_TRUE(contains(hash_set, 42));
    EXPECT_FALSE(contains(hash_set, 7));
    EXPECT_FALSE(contains(hash_set, -42));
}

TEST(Contains, SetIntEmpty) {
    std::unordered_set<int> hash_set;
    EXPECT_FALSE(contains(hash_set, 0));
}

TEST(Contains, SetString) {
    std::unordered_set<std::string> hash_set = {"", "ACGT"};
    EXPECT_TRUE(contains(hash_set, std::string("")));
    EXPECT_TRUE(contains(hash_set, std::string("ACGT")));
    // matches are exact: lowercase, prefixes, and extensions are all absent
    EXPECT_FALSE(contains(hash_set, std::string("acgt")));
    EXPECT_FALSE(contains(hash_set, std::string("ACG")));
    EXPECT_FALSE(contains(hash_set, std::string("ACGTA")));
}

// The unordered_set<idx4> instantiation below is the one calc_prec_recall_aln itself uses.

TEST(Contains, SetPresent) {
    std::unordered_set<idx4> wave = {idx4(0,0,0,0), idx4(1,2,3,4), idx4(5,0,0,7)};
    EXPECT_TRUE(contains(wave, idx4(0,0,0,0)));
    EXPECT_TRUE(contains(wave, idx4(1,2,3,4)));
    EXPECT_TRUE(contains(wave, idx4(5,0,0,7)));
}

TEST(Contains, SetAbsent) {
    std::unordered_set<idx4> wave = {idx4(1,2,3,4)};
    EXPECT_FALSE(contains(wave, idx4(1,2,3,5)));
    EXPECT_FALSE(contains(wave, idx4(0,0,0,0)));
}

TEST(Contains, SetIdx4DistinguishesEachField) {
    // the hash XORs one term per field, so a lookup must still be decided by idx4::operator==
    std::unordered_set<idx4> wave = {idx4(1,0,0,0)};
    EXPECT_TRUE(contains(wave, idx4(1,0,0,0)));
    EXPECT_FALSE(contains(wave, idx4(0,1,0,0)));
    EXPECT_FALSE(contains(wave, idx4(0,0,1,0)));
    EXPECT_FALSE(contains(wave, idx4(0,0,0,1)));
}

/* calc_ng50 **************************************************************************************/

TEST(TestNG50, TestNG50Calc) {
    std::vector<int> phase_blocks1 = {25, 50};
    EXPECT_EQ(50, calc_ng50(phase_blocks1, 100));

    std::vector<int> phase_blocks2 = {1, 25, 25, 50};
    EXPECT_EQ(25, calc_ng50(phase_blocks2, 101));

    std::vector<int> phase_blocks3 = {90};
    EXPECT_EQ(90, calc_ng50(phase_blocks3, 100));

    std::vector<int> phase_blocks4 = {1, 10, 20};
    EXPECT_EQ(0, calc_ng50(phase_blocks4, 100));
}

TEST(CalcNg50, Empty) {
    // no blocks means nothing is phased, so the loop never runs and the fallback 0 is returned
    std::vector<int> phase_blocks;
    EXPECT_EQ(0, calc_ng50(phase_blocks, 100));
}

TEST(CalcNg50, TotalZero) {
    // with a zero genome size the 50% threshold is 0, met by the largest block on the first step
    std::vector<int> phase_blocks = {10, 20};
    EXPECT_EQ(20, calc_ng50(phase_blocks, 0));
}

TEST(CalcNg50, EmptyAndZero) {
    // the threshold is vacuously met, but with no blocks there is nothing to return
    std::vector<int> phase_blocks;
    EXPECT_EQ(0, calc_ng50(phase_blocks, 0));
}

TEST(CalcNg50, SingleBelow) {
    // 40 < 50, so the single block never reaches half the genome
    std::vector<int> phase_blocks = {40};
    EXPECT_EQ(0, calc_ng50(phase_blocks, 100));
}

TEST(CalcNg50, ExactEvenTies) {
    // four equal blocks: the running sum first reaches 50 on the second block, which is also 25
    std::vector<int> phase_blocks = {25, 25, 25, 25};
    EXPECT_EQ(25, calc_ng50(phase_blocks, 100));
}

TEST(CalcNg50, UnsortedInput) {
    // sorted descending to {90, 50, 10}; 90 < 100 but 90+50 >= 100, so NG50 is 50
    std::vector<int> phase_blocks = {10, 90, 50};
    EXPECT_EQ(50, calc_ng50(phase_blocks, 200));

    // every permutation of the same multiset gives the same answer
    std::vector<int> permuted = {50, 10, 90};
    EXPECT_EQ(50, calc_ng50(permuted, 200));
    std::vector<int> descending = {90, 50, 10};
    EXPECT_EQ(50, calc_ng50(descending, 200));
}

TEST(CalcNg50, InputUnmodified) {
    // phase_blocks is taken by value, so the caller's ordering survives the internal sort
    std::vector<int> phase_blocks = {10, 90, 50};
    calc_ng50(phase_blocks, 200);
    EXPECT_EQ(std::vector<int>({10, 90, 50}), phase_blocks);
}

/* generate_str ***********************************************************************************/

TEST(GenerateStr, NoVariants) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGTACGT");
    auto vars = make_ctgVariants("chr1", {});

    // with no variants the result is exactly the half-open reference slice [2, 6)
    EXPECT_EQ("GTAC", generate_str(ref, vars, "chr1", 0, 0, 2, 6));
}

TEST(GenerateStr, Snp) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "T"}});

    // "AC" + alt "T" + ref[3,8) "TACGT"
    EXPECT_EQ("ACTTACGT", generate_str(ref, vars, "chr1", 0, 1, 0, 8));
}

TEST(GenerateStr, Ins) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{2, 0, TYPE_INS, "", "TT"}});

    // an insertion consumes no reference: "AC" + alt "TT" + ref[2,8) "GTACGT"
    EXPECT_EQ("ACTTGTACGT", generate_str(ref, vars, "chr1", 0, 1, 0, 8));
}

TEST(GenerateStr, Del) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{2, 2, TYPE_DEL, "GT", ""}});

    // a deletion emits nothing and skips its two reference bases: "AC" + ref[4,8) "ACGT"
    EXPECT_EQ("ACACGT", generate_str(ref, vars, "chr1", 0, 1, 0, 8));
}

TEST(GenerateStr, Cpx) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{2, 2, TYPE_CPX, "GT", "AAA"}});

    // a complex variant emits its alt and skips its whole ref allele: "AC" + "AAA" + ref[4,8)
    EXPECT_EQ("ACAAAACGT", generate_str(ref, vars, "chr1", 0, 1, 0, 8));
}

TEST(GenerateStr, MinQualFiltersSub) {
    GlobalsGuard guard;
    g.max_qual = 60;
    g.min_qual = 20;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "T", GT_REF_ALT1, 10}});

    // the filtered substitution advances neither the sequence nor ref_pos, so its ref base is kept
    EXPECT_EQ("ACGTACGT", generate_str(ref, vars, "chr1", 0, 1, 0, 8, g.min_qual));
}

TEST(GenerateStr, MinQualFiltersDel) {
    GlobalsGuard guard;
    g.max_qual = 60;
    g.min_qual = 20;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{2, 2, TYPE_DEL, "GT", "", GT_REF_ALT1, 10}});

    // a filtered deletion does not skip its reference bases, so "GT" is still emitted
    EXPECT_EQ("ACGTACGT", generate_str(ref, vars, "chr1", 0, 1, 0, 8, g.min_qual));
}

TEST(GenerateStr, MinQualKeepsAtThreshold) {
    GlobalsGuard guard;
    g.max_qual = 60;
    g.min_qual = 20;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "T", GT_REF_ALT1, 20}});

    // the comparison is >=, so a variant exactly at the threshold is applied
    EXPECT_EQ("ACTTACGT", generate_str(ref, vars, "chr1", 0, 1, 0, 8, g.min_qual));
}

TEST(GenerateStr, PrefixSkip) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1",
            {{1, 1, TYPE_SUB, "C", "T"}, {5, 1, TYPE_SUB, "C", "A"}});

    // the variant at position 1 is skipped as < beg_pos: ref[3,5) "TA" + alt "A" + ref[6,8) "GT"
    EXPECT_EQ("TAAGT", generate_str(ref, vars, "chr1", 0, 2, 3, 8));
}

TEST(GenerateStr, BoundaryEnd) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {{5, 1, TYPE_SUB, "C", "A"}});

    // end_pos is exclusive, so a variant at position 5 is outside [0, 5) and never applied
    EXPECT_EQ("ACGTA", generate_str(ref, vars, "chr1", 0, 1, 0, 5));
}

TEST(GenerateStr, MissingContigErrors) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr2", {});

    EXPECT_EXIT(generate_str(ref, vars, "chr2", 0, 0, 0, 4),
            testing::ExitedWithCode(1), "not in reference FASTA");
}

TEST(GenerateStr, PosPastContigEndErrors) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {});

    // the contig is present, so the position is named rather than blamed on the contig (#167)
    EXPECT_EXIT(generate_str(ref, vars, "chr1", 0, 0, 10, 12),
            testing::ExitedWithCode(1), "Position 10 out of range on contig 'chr1' of length 8");
}

TEST(GenerateStr, NegativePosErrors) {
    GlobalsGuard guard;
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto vars = make_ctgVariants("chr1", {});

    // a negative position would otherwise convert to a huge size_t inside substr() (#167)
    EXPECT_EXIT(generate_str(ref, vars, "chr1", 0, 0, -1, 2),
            testing::ExitedWithCode(1), "Position -1 out of range on contig 'chr1' of length 8");
}

/* wf_ed ******************************************************************************************/

TEST(WfEd, Identical) {
    int score = -1;
    wf_ed("ACGT", "ACGT", score);
    EXPECT_EQ(0, score);
}

TEST(WfEd, SingleSub) {
    // one mismatched base at index 1
    int score = -1;
    wf_ed("ACGT", "AGGT", score);
    EXPECT_EQ(1, score);
}

TEST(WfEd, SingleIns) {
    // the query carries one extra G
    int score = -1;
    wf_ed("ACGGT", "ACGT", score);
    EXPECT_EQ(1, score);
}

TEST(WfEd, SingleDel) {
    // the query is missing the truth's G
    int score = -1;
    wf_ed("ACT", "ACGT", score);
    EXPECT_EQ(1, score);
}

TEST(WfEd, EmptyQuery) {
    // an empty query costs one insertion per truth base
    int score = -1;
    wf_ed("", "ACGT", score);
    EXPECT_EQ(4, score);
}

TEST(WfEd, EmptyTruth) {
    // an empty truth costs one deletion per query base
    int score = -1;
    wf_ed("ACGT", "", score);
    EXPECT_EQ(4, score);
}

TEST(WfEd, BothEmpty) {
    int score = -1;
    wf_ed("", "", score);
    EXPECT_EQ(0, score);
}

TEST(WfEd, TwoSubs) {
    // mismatches at indices 1 (C/G) and 6 (G/C)
    int score = -1;
    wf_ed("ACGTACGT", "AGGTACCT", score);
    EXPECT_EQ(2, score);
}

TEST(WfEd, AllDifferent) {
    // equal lengths with no shared subsequence, so every base must be substituted
    int score = -1;
    wf_ed("AAAA", "CCCC", score);
    EXPECT_EQ(4, score);
}

TEST(WfEd, ShiftedRepeat) {
    // both are 8bp and every column mismatches, so no 1-edit alignment exists; deleting the
    // trailing T and inserting a leading T costs 2
    int score = -1;
    wf_ed("ACGTACGT", "TACGTACG", score);
    EXPECT_EQ(2, score);
}

TEST(WfEd, Symmetry) {
    // unit-cost edit distance is a metric, so swapping the arguments cannot change it
    const std::vector< std::pair<std::string, std::string> > pairs = {
        {"ACGT", "ACGT"}, {"ACGT", "AGGT"}, {"ACGGT", "ACGT"}, {"ACT", "ACGT"},
        {"", "ACGT"}, {"AAAA", "CCCC"}, {"ACGTACGT", "TACGTACG"}, {"GATTACA", "GATACA"}};
    for (const auto & [a, b] : pairs) {
        int forward = -1;
        int backward = -1;
        wf_ed(a, b, forward);
        wf_ed(b, a, backward);
        EXPECT_EQ(forward, backward) << "'" << a << "' vs '" << b << "'";
    }
}

/* wf_swg_align ***********************************************************************************/

// Every case that reaches the wavefront proper keeps both strings at least one character long.
// The wavefront itself indexes a mat_len of query_len+truth_len-1 and would read out of range on
// an empty input, so it is only the early-exit guard added for issue #65 that makes an empty
// string safe; staying inside the supported domain is deliberate rather than an oversight, and
// the three guarded shapes are asserted separately below.

TEST(WfSwgAlign, Identical) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    int score = -1;
    wf_swg_align("ACGT", "ACGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(0, score);
}

TEST(WfSwgAlign, SingleSub) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // one substitution costs x=4; replacing it with a deletion plus an insertion costs
    // 2*(o+e) = 16, so the substitution wins
    int score = -1;
    wf_swg_align("ACGT", "AGGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(g.sub, score);
}

TEST(WfSwgAlign, GapLen1) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // the query is missing one G, and "ACGT" embeds in "ACGGT" with no mismatches, so the whole
    // cost is a single length-1 gap: o + e = 8
    int score = -1;
    wf_swg_align("ACGT", "ACGGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(g.open + g.extend, score);
}

TEST(WfSwgAlign, GapLen1Reversed) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // the same length-1 gap, with the extra base on the query side
    int score = -1;
    wf_swg_align("ACGGT", "ACGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(g.open + g.extend, score);
}

TEST(WfSwgAlign, GapLen2Affine) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // "ACGT" embeds in "ACGGGT" by removing the two contiguous middle Gs, so one length-2 gap
    // with no mismatches suffices. Gap-affine charges o + 2e = 10 for that, NOT 2*(o+e) = 16;
    // this is the assertion a broken affine term fails.
    int score = -1;
    wf_swg_align("ACGT", "ACGGGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(g.open + 2*g.extend, score);
    EXPECT_NE(2*(g.open + g.extend), score);
}

TEST(WfSwgAlign, TwoSeparateGaps) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // "ACGTACGT" embeds in "ACGGTACGGT" only by removing two non-adjacent Gs, so two length-1
    // gaps are needed: 2*(o+e) = 16. No contiguous 2-base removal from the truth leaves fewer
    // than 3 mismatches, so the single-gap alternative costs at least (o+2e) + 3x = 22.
    int score = -1;
    wf_swg_align("ACGTACGT", "ACGGTACGGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(2*(g.open + g.extend), score);
}

TEST(WfSwgAlign, SubVsGapChoiceSubWins) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // x=4 is cheaper than the 2*(o+e)=16 deletion-plus-insertion detour
    int score = -1;
    wf_swg_align("ACGT", "AGGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(4, score);
}

TEST(WfSwgAlign, SubVsGapChoiceGapWins) {
    GlobalsGuard guard;
    set_penalties(20, 2, 1);

    // same strings, but now x=20 exceeds the 2*(o+e)=6 detour, so the gap pair is chosen
    int score = -1;
    wf_swg_align("ACGT", "AGGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(6, score);
}

TEST(WfSwgAlign, BothEmpty) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // two empty sequences are already aligned, so no gap is opened
    int score = -1;
    wf_swg_align("", "", score, g.sub, g.open, g.extend);
    EXPECT_EQ(0, score);
}

TEST(WfSwgAlign, EmptyQuery) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // the whole truth is one gap: o + e*4 = 14
    int score = -1;
    wf_swg_align("", "ACGT", score, g.sub, g.open, g.extend);
    EXPECT_EQ(g.open + 4*g.extend, score);
}

TEST(WfSwgAlign, EmptyTruth) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // symmetrically, the whole query is one gap: o + e*4 = 14
    int score = -1;
    wf_swg_align("ACGT", "", score, g.sub, g.open, g.extend);
    EXPECT_EQ(g.open + 4*g.extend, score);
}

TEST(WfSwgAlign, PenaltySensitivity) {
    GlobalsGuard guard;

    // for a single mismatch the score is exactly x, as long as x stays below the 2*(o+e) detour
    for (int sub : {2, 4, 8, 15}) {
        set_penalties(sub, 6, 2);
        int score = -1;
        wf_swg_align("ACGT", "AGGT", score, g.sub, g.open, g.extend);
        EXPECT_EQ(sub, score) << "x = " << sub;
    }
}

/* wf_swg_max_reach *******************************************************************************/

// The main diagonal guard only bites when main_diag_off is small, so unless a test is exercising
// it, main_diag_start is set far beyond either sequence.
static const int kNoMainDiagBlock = 1000;

/**
 * @brief Allocates a fresh offsets buffer for one wf_swg_max_reach call over the given sequences.
 *
 * wf_swg_max_reach requires its caller to allocate the buffer and reads stale entries, so every
 * call needs its own. Sizing is delegated to alloc_reach_offs and the penalties come from `g`.
 */
static std::vector<int> reach_offs(const std::string & query, const std::string & truth) {
    return alloc_reach_offs(int(query.size()), int(truth.size()), g.sub, g.open, g.extend);
}

TEST(WfSwgMaxReach, FullMatch) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGT";
    const std::string truth = "ACGT";
    std::vector<int> offs = reach_offs(query, truth);

    // identical strings extend to the final column at score 0, so the last truth index is reached
    EXPECT_EQ(3, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, 0,
                g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, PartialBudget) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTAAA";
    const std::string truth = "ACGTCCC";

    // the shared prefix "ACGT" extends for free, then every remaining column mismatches, so a
    // zero budget stops at truth index 3
    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_EQ(3, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, 0,
                g.sub, g.open, g.extend));

    // two substitutions (2x = 8) buy two more bases of the mismatched tail
    std::vector<int> offs2 = reach_offs(query, truth);
    EXPECT_EQ(5, wf_swg_max_reach(query, truth, offs2, 0, kNoMainDiagBlock, 8,
                g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, ScoreAllowsSub) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTAAA";
    const std::string truth = "ACGTCCC";

    // each substitution costs x=4 and buys exactly one more base of reach past the shared prefix;
    // budgets between multiples of x buy nothing, since no gap is affordable below o+e=8 either
    const std::vector< std::pair<int,int> > budget_to_reach =
        {{0, 3}, {3, 3}, {4, 4}, {7, 4}, {8, 5}};
    for (const auto & [max_score, expected] : budget_to_reach) {
        std::vector<int> offs = reach_offs(query, truth);
        EXPECT_EQ(expected, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, max_score,
                    g.sub, g.open, g.extend)) << "max_score = " << max_score;
    }
}

TEST(WfSwgMaxReach, PrefersBetterDiagonalWhenQueryExhausted) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTAAA";
    const std::string truth = "ACGTCCC";

    // a budget of 12 affords three substitutions along the main diagonal (3x = 12), reaching the
    // final truth index 6. It also affords an insertion run (o+e, then e per extension: 8, 10, 12)
    // that consumes all three remaining query bases without advancing along the truth, reaching
    // only index 3. The extend loop scans diagonals in ascending index order and so meets the
    // insertion-heavy diagonal first, but the reported reach must be the best over all diagonals.
    std::vector<int> offs11 = reach_offs(query, truth);
    EXPECT_EQ(5, wf_swg_max_reach(query, truth, offs11, 0, kNoMainDiagBlock, 11,
                g.sub, g.open, g.extend));
    std::vector<int> offs12 = reach_offs(query, truth);
    EXPECT_EQ(6, wf_swg_max_reach(query, truth, offs12, 0, kNoMainDiagBlock, 12,
                g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, MonotonicInScoreBudget) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTAAA";
    const std::string truth = "ACGTCCC";

    // a larger budget can only afford a superset of the alignments a smaller one affords, so reach
    // must never decrease as max_score grows
    int prev = 0;
    for (int max_score = 0; max_score <= 24; max_score++) {
        std::vector<int> offs = reach_offs(query, truth);
        int reach = wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, max_score,
                g.sub, g.open, g.extend);
        EXPECT_LE(prev, reach) << "reach decreased at max_score = " << max_score;
        prev = reach;
    }
}

TEST(WfSwgMaxReach, MonotonicInScoreBudgetReverse) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTAAA";
    const std::string truth = "ACGTCCC";

    // the reverse variant leaves a gap for o rather than opening one forwards for o+e, and must
    // satisfy the same monotonicity
    int prev = 0;
    for (int max_score = 0; max_score <= 24; max_score++) {
        std::vector<int> offs = reach_offs(query, truth);
        int reach = wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, max_score,
                g.sub, g.open, g.extend, false /* print */, true /* reverse */);
        EXPECT_LE(prev, reach) << "reach decreased at max_score = " << max_score;
        prev = reach;
    }
}

TEST(WfSwgMaxReach, MonotonicWhenDeletionsAdvanceFarthest) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "AAAA";
    const std::string truth = "ACGTACGTA";

    // the query is shorter than the truth, so advancing past index 0 needs deletions rather than
    // substitutions: query[0] matches truth[0] for free, then deleting "CGT" costs o + 3e = 12 and
    // query[1] == 'A' matches truth[4], for a reach of 4
    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_EQ(4, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, 12,
                g.sub, g.open, g.extend));

    int prev = 0;
    for (int max_score = 0; max_score <= 24; max_score++) {
        std::vector<int> sweep_offs = reach_offs(query, truth);
        int reach = wf_swg_max_reach(query, truth, sweep_offs, 0, kNoMainDiagBlock, max_score,
                g.sub, g.open, g.extend);
        EXPECT_LE(prev, reach) << "reach decreased at max_score = " << max_score;
        prev = reach;
    }
}

TEST(WfSwgMaxReach, GapExtensionAffine) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTACGT";
    const std::string truth = "ACGTGGACGT";

    // the truth has two extra bases, so reaching its end needs a length-2 deletion: o + 2e = 10
    // under gap-affine scoring, versus 2*(o+e) = 16 if extension were charged as a second open.
    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_EQ(9, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, g.open + 2*g.extend,
                g.sub, g.open, g.extend));

    // one point short of o+2e the gap cannot be extended, and the best remaining path is two
    // substitutions along the main diagonal, reaching truth index 5
    std::vector<int> offs2 = reach_offs(query, truth);
    EXPECT_EQ(5, wf_swg_max_reach(query, truth, offs2, 0, kNoMainDiagBlock,
                g.open + 2*g.extend - 1, g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, ReverseFlag) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTACGT";
    const std::string truth = "ACGTGACGT";

    // forwards a gap is entered for o+e=8, so a budget of 8 reaches the last truth index while a
    // budget of 7 is stuck at the reach of a single substitution
    std::vector<int> fwd8 = reach_offs(query, truth);
    EXPECT_EQ(8, wf_swg_max_reach(query, truth, fwd8, 0, kNoMainDiagBlock, 8,
                g.sub, g.open, g.extend, /* print = */ false, /* reverse = */ false));
    std::vector<int> fwd7 = reach_offs(query, truth);
    EXPECT_EQ(4, wf_swg_max_reach(query, truth, fwd7, 0, kNoMainDiagBlock, 7,
                g.sub, g.open, g.extend, /* print = */ false, /* reverse = */ false));

    // in reverse the gap-open penalty is charged on leaving the gap instead of entering it, so a
    // budget of 7 already pays for three extension steps and reaches truth index 6
    std::vector<int> rev8 = reach_offs(query, truth);
    EXPECT_EQ(8, wf_swg_max_reach(query, truth, rev8, 0, kNoMainDiagBlock, 8,
                g.sub, g.open, g.extend, /* print = */ false, /* reverse = */ true));
    std::vector<int> rev7 = reach_offs(query, truth);
    EXPECT_EQ(6, wf_swg_max_reach(query, truth, rev7, 0, kNoMainDiagBlock, 7,
                g.sub, g.open, g.extend, /* print = */ false, /* reverse = */ true));
}

TEST(WfSwgMaxReach, MainDiagBlock) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTACGT";
    const std::string truth = "ACGTACGT";

    // main_diag_off = main_diag_start - main_diag = 3, and the guard stops diagonal extension
    // once off+1 reaches it, so the match run ends at query offset 2
    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_EQ(2, wf_swg_max_reach(query, truth, offs, 0, 3, 0, g.sub, g.open, g.extend));

    // without the block the identical strings extend all the way to the final truth index
    std::vector<int> offs2 = reach_offs(query, truth);
    EXPECT_EQ(7, wf_swg_max_reach(query, truth, offs2, 0, kNoMainDiagBlock, 0,
                g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, MainDiagBlockOnlyBlocksMainDiag) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTACGT";
    const std::string truth = "ACGTACGT";

    // same main_diag_off = 3, but declared on diagonal 1; the wavefront starts on diagonal 0 and
    // so is unaffected
    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_EQ(7, wf_swg_max_reach(query, truth, offs, 1, 4, 0, g.sub, g.open, g.extend));
}

// An empty sequence reaches no truth base, so the guarded return of 0 is a sentinel rather than a
// real index; it is indistinguishable from genuinely reaching truth index 0, and matches the value
// max_reach already falls through with when no diagonal holds a valid offset.
//
// The empty-query cases are the memory-safety regression tests: without the guard the seed write at
// diagonal query_len-1 = -1 lands before the caller's buffer. Empty truth alone stays in bounds but
// misplaces the seed into the next score row, and the value that leaks out depends on how many
// score rows are visited before s2 wraps -- at x=4, o=6, e=2 it was 0 up to max_score 12 and -1
// from 13 on. So that case pins a value where there was previously no consistent one.

TEST(WfSwgMaxReach, EmptyQuery) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "";
    const std::string truth = "ACGT";

    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_EQ(0, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, 1000,
                g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, EmptyTruth) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGT";
    const std::string truth = "";

    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_EQ(0, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, 1000,
                g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, BothEmpty) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "";
    const std::string truth = "";

    // the buffer is empty here, so any indexing at all is out of bounds
    std::vector<int> offs = reach_offs(query, truth);
    EXPECT_TRUE(offs.empty());
    EXPECT_EQ(0, wf_swg_max_reach(query, truth, offs, 0, kNoMainDiagBlock, 1000,
                g.sub, g.open, g.extend));
}

TEST(WfSwgMaxReach, EmptyReverse) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);

    // the guard precedes the forward/reverse split, so it holds in reverse too
    std::vector<int> offs = reach_offs("", "ACGT");
    EXPECT_EQ(0, wf_swg_max_reach("", "ACGT", offs, 0, kNoMainDiagBlock, 1000,
                g.sub, g.open, g.extend, /* print = */ false, /* reverse = */ true));
}

/* idx4 *******************************************************************************************/

/** @brief Returns idx4 values covering every single-field variation used by the ordering tests. */
static std::vector<idx4> ordering_samples() {
    return {idx4(0,0,0,0), idx4(0,0,0,1), idx4(0,0,1,0), idx4(0,0,1,1),
            idx4(0,1,0,0), idx4(0,1,2,3), idx4(1,0,0,0), idx4(1,0,0,5),
            idx4(1,2,3,4), idx4(2,0,0,0), idx4(-1,0,0,0), idx4(0,-1,0,0),
            idx4(0,0,-1,0), idx4(0,0,0,-1)};
}

TEST(Idx4, DefaultCtor) {
    idx4 a;
    EXPECT_EQ(0, a.qni);
    EXPECT_EQ(0, a.tni);
    EXPECT_EQ(0, a.qi);
    EXPECT_EQ(0, a.ti);
}

TEST(Idx4, EqNeq) {
    const idx4 a(1,2,3,4);
    EXPECT_TRUE(a == idx4(1,2,3,4));
    EXPECT_FALSE(a != idx4(1,2,3,4));

    // a difference in any single field breaks equality
    EXPECT_TRUE(a != idx4(9,2,3,4));
    EXPECT_TRUE(a != idx4(1,9,3,4));
    EXPECT_TRUE(a != idx4(1,2,9,4));
    EXPECT_TRUE(a != idx4(1,2,3,9));
    EXPECT_FALSE(a == idx4(9,2,3,4));
    EXPECT_FALSE(a == idx4(1,9,3,4));
    EXPECT_FALSE(a == idx4(1,2,9,4));
    EXPECT_FALSE(a == idx4(1,2,3,9));
}

TEST(Idx4, CopyAndAssign) {
    const idx4 a(1,2,3,4);
    idx4 copied(a);
    EXPECT_TRUE(copied == a);

    idx4 assigned;
    assigned = a;
    EXPECT_TRUE(assigned == a);

    // self-assignment is a no-op rather than a corruption (indirected to avoid -Wself-assign)
    idx4 * self = &assigned;
    assigned = *self;
    EXPECT_TRUE(assigned == a);
}

TEST(Idx4, HashEqualForEqualKeys) {
    const std::hash<idx4> hasher;
    EXPECT_EQ(hasher(idx4(1,2,3,4)), hasher(idx4(1,2,3,4)));
    EXPECT_EQ(hasher(idx4()), hasher(idx4(0,0,0,0)));
}

TEST(Idx4, HashUsableInSet) {
    std::unordered_set<idx4> visited;
    visited.insert(idx4(1,2,3,4));
    visited.insert(idx4(1,2,3,4));
    EXPECT_EQ(size_t(1), visited.size());

    visited.insert(idx4(4,3,2,1));
    EXPECT_EQ(size_t(2), visited.size());
    EXPECT_EQ(size_t(1), visited.count(idx4(1,2,3,4)));
    EXPECT_EQ(size_t(0), visited.count(idx4(1,2,3,5)));
}

TEST(Idx4, HashUsableInMap) {
    // a repeated key overwrites rather than adding an entry, as calc_prec_recall_aln relies on
    std::unordered_map<idx4, idx4> ptrs;
    ptrs[idx4(7,8,9,10)] = idx4(1,2,3,3);
    ptrs[idx4(7,8,9,10)] = idx4(1,2,2,3);
    EXPECT_EQ(size_t(1), ptrs.size());
    EXPECT_TRUE(idx4(1,2,2,3) == ptrs.at(idx4(7,8,9,10)));

    // every single-field variation is a distinct key, none of which collides with idx4(7,8,9,10)
    const std::vector<idx4> samples = ordering_samples();
    for (const idx4 & sample : samples) ptrs[sample] = idx4();
    EXPECT_EQ(samples.size() + 1, ptrs.size());
}

TEST(Idx4, LessLexicographic) {
    // qni dominates, then tni, then qi, then ti
    EXPECT_TRUE(idx4(0,9,9,9) < idx4(1,0,0,0));
    EXPECT_TRUE(idx4(1,0,9,9) < idx4(1,1,0,0));
    EXPECT_TRUE(idx4(1,1,0,9) < idx4(1,1,1,0));
    EXPECT_TRUE(idx4(1,1,1,0) < idx4(1,1,1,1));
    EXPECT_FALSE(idx4(1,0,0,0) < idx4(0,9,9,9));
    EXPECT_FALSE(idx4(1,1,1,1) < idx4(1,1,1,1));
}

TEST(Idx4, LessStrictWeakOrdering) {
    // #66 made operator< a valid strict weak ordering, so these properties are enforced here
    const std::vector<idx4> samples = ordering_samples();

    // irreflexivity: no value precedes itself
    for (const idx4 & a : samples) EXPECT_FALSE(a < a);

    // asymmetry: at most one of a<b and b<a holds
    for (const idx4 & a : samples) {
        for (const idx4 & b : samples) {
            if (a < b) EXPECT_FALSE(b < a);
        }
    }

    // transitivity of the order, and of the induced equivalence
    for (const idx4 & a : samples) {
        for (const idx4 & b : samples) {
            for (const idx4 & c : samples) {
                if (a < b && b < c) EXPECT_TRUE(a < c);
                const bool ab_equiv = !(a < b) && !(b < a);
                const bool bc_equiv = !(b < c) && !(c < b);
                if (ab_equiv && bc_equiv) {
                    EXPECT_TRUE(!(a < c) && !(c < a));
                }
            }
        }
    }

    // the induced equivalence is exactly operator==, so a total order on distinct values
    for (const idx4 & a : samples) {
        for (const idx4 & b : samples) {
            const bool equiv = !(a < b) && !(b < a);
            EXPECT_EQ(a == b, equiv);
        }
    }
}

/* Graph fixtures *********************************************************************************/

/**
 * @struct GraphFixture
 * @brief A reference, the query and truth variant containers over it, and the graph they build.
 *
 * Every member is held so that a test can assert on the graph and then read the labels the
 * evaluation writes back into the variant containers.
 */
struct GraphFixture {
    std::shared_ptr<fastaData> ref;        ///< Single-contig reference named "chr1"
    std::shared_ptr<ctgVariants> qvars;    ///< Query variants
    std::shared_ptr<ctgVariants> tvars;    ///< Truth variants
    std::shared_ptr<ctgSuperclusters> sc;  ///< Supercluster holding both containers
    std::shared_ptr<Graph> graph;          ///< Graph over supercluster 0
};

/**
 * @brief Builds a one-contig, one-supercluster graph fixture over the given reference and variants.
 *
 * Every variant is forced into supercluster 0, because the Graph constructor takes its variant
 * range from lower_bound/upper_bound over `superclusters` and the var_desc default of -1 would
 * silently leave the fixture with no variants at all.
 */
static GraphFixture build_fixture(const std::string & ref_seq, std::vector<var_desc> qv,
        std::vector<var_desc> tv, int truth_hap = HAP1) {
    for (var_desc & v : qv) v.supercluster = 0;
    for (var_desc & v : tv) v.supercluster = 0;
    GraphFixture f;
    f.ref = make_fasta("chr1", ref_seq);
    f.qvars = make_ctgVariants("chr1", qv);
    f.tvars = make_ctgVariants("chr1", tv);
    f.sc = make_ctgSuperclusters(f.qvars, f.tvars);
    f.graph = make_graph(f.sc, f.ref, "chr1", truth_hap);
    return f;
}

/** @brief Aligns and labels one graph exactly as evaluate_variants does, returning the score. */
static int align_and_label(std::shared_ptr<Graph> graph, int truth_hap) {
    std::unordered_map<idx4, idx4> ptrs;
    int score = calc_prec_recall_aln(graph, ptrs, false);
    calc_prec_recall(graph, ptrs, truth_hap, false);
    return score;
}

/** @brief Returns a length-len reference cycling through ACGT, for fixtures needing a long ref. */
static std::string periodic_ref(int len) {
    std::string seq;
    for (int i = 0; i < len; i++) seq += "ACGT"[i % 4];
    return seq;
}

/* Graph ******************************************************************************************/

// The fixtures below always carry at least one variant on one callset. A supercluster empty on
// both is not a defined input: get_min_ref_pos() returns int::max() - 1 for it, and the
// constructor takes a reference substring at that offset.

TEST(GraphCtor, QuerySnpNodeLayout) {
    // ref "ACGTACGT" with one query SNP at pos 2 (G->T) and no truth variants. The span is
    // [pos-1, pos+rlen+1] = [1, 4], so nodes are cut at reference offsets 1, 2, 3 and 5 (the
    // trailing node runs one past ref_end), giving offsets 0, 1, 2 and 4 relative to ref_beg.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 60, 0, 0}}, {});

    ASSERT_EQ(4, f.graph->qnodes);
    EXPECT_EQ(std::vector<std::string>({"_C", "_T", "_G", "_TA"}), f.graph->qseqs);
    EXPECT_EQ(std::vector<int>({0, 1, 1, 2}), f.graph->qbegs);
    EXPECT_EQ(std::vector<int>({1, 2, 2, 4}), f.graph->qends);
    EXPECT_EQ(std::vector<int>({TYPE_REF, TYPE_SUB, TYPE_REF, TYPE_REF}), f.graph->qtypes);
    EXPECT_EQ(std::vector<int>({-1, 0, -1, -1}), f.graph->qidxs);

    // no truth variants, so the truth side is a single reference node spanning the whole window
    ASSERT_EQ(1, f.graph->tnodes);
    EXPECT_EQ("_CGTA", f.graph->tseqs[0]);
    EXPECT_EQ(-1, f.graph->tskips[0]);
}

TEST(GraphCtor, ContigStartRefSpanClamped) {
    // A variant at position 0 has no base to its left, so ref_beg clamps to 0 instead of -1 and
    // the constructor no longer slices the reference at a negative offset (#167).
    GraphFixture f = build_fixture("ACGTACGT",
            {{0, 1, TYPE_SUB, "A", "G", GT_ALT1_REF, 60, 0, 0}}, {});

    EXPECT_EQ("ACG", f.graph->ref);
    EXPECT_EQ(std::vector<int>({0, 0, 1}), f.graph->qbegs);
    EXPECT_EQ(std::vector<int>({1, 1, 3}), f.graph->qends);
}

TEST(GraphCtor, ContigStartHasNoLeftFlankNode) {
    // Away from the contig start the first node is the one-base left flank (QuerySnpNodeLayout).
    // At position 0 that node cannot exist, so the graph opens directly on the variant and its
    // parallel reference allele, which no edge reaches from node 0.
    GraphFixture f = build_fixture("ACGTACGT",
            {{0, 1, TYPE_SUB, "A", "G", GT_ALT1_REF, 60, 0, 0}}, {});

    ASSERT_EQ(3, f.graph->qnodes);
    EXPECT_EQ(std::vector<std::string>({"_G", "_A", "_CG"}), f.graph->qseqs);
    EXPECT_EQ(std::vector<int>({TYPE_SUB, TYPE_REF, TYPE_REF}), f.graph->qtypes);
    EXPECT_TRUE(f.graph->qprevs[0].empty());
    EXPECT_TRUE(f.graph->qprevs[1].empty()) << "the reference allele is unreachable from node 0";
}

TEST(GraphCtor, QueryVariantHasParallelRefAllele) {
    // Adding a variant node does not advance ref_pos; the variant's end is pushed onto
    // qnode_ends and popped later, so the reference allele reappears as node 2 spanning the same
    // [1, 2) the SNP does. The query graph therefore offers both alleles with no bypass machinery.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 60, 0, 0}}, {});

    ASSERT_EQ(4, f.graph->qnodes);
    EXPECT_EQ(f.graph->qbegs[1], f.graph->qbegs[2]);
    EXPECT_EQ(f.graph->qends[1], f.graph->qends[2]);
    EXPECT_EQ(TYPE_SUB, f.graph->qtypes[1]);
    EXPECT_EQ(TYPE_REF, f.graph->qtypes[2]);
    EXPECT_EQ("_G", f.graph->qseqs[2]) << "node 2 should carry the reference allele";
}

TEST(GraphCtor, RefSpanFromVariantBounds) {
    // ref_beg = min_pos - 1 and ref_end = max(pos + rlen) + 1, and the stored slice runs to
    // ref_end + 1 exclusive: one base of left flank and two bases past the variant's end.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 60, 0, 0}}, {});

    EXPECT_EQ("CGTA", f.graph->ref);
    EXPECT_EQ(std::string("ACGTACGT").substr(1, 4), f.graph->ref);
}

TEST(GraphCtor, QueryPointersWireParallelAlleles) {
    // node 0 [0,1) precedes both alleles at [1,2), and both lead into node 3 at [2,4).
    // Successor lists are built by an ascending scan, so they are in ascending node order.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 60, 0, 0}}, {});

    ASSERT_EQ(4, f.graph->qnodes);
    EXPECT_EQ(std::vector<int>({1, 2}), f.graph->qnexts[0]);
    EXPECT_EQ(std::vector<int>({3}), f.graph->qnexts[1]);
    EXPECT_EQ(std::vector<int>({3}), f.graph->qnexts[2]);
    EXPECT_TRUE(f.graph->qnexts[3].empty());
    EXPECT_TRUE(f.graph->qprevs[0].empty());
    EXPECT_EQ(std::vector<int>({0}), f.graph->qprevs[1]);
    EXPECT_EQ(std::vector<int>({0}), f.graph->qprevs[2]);
    EXPECT_EQ(std::vector<int>({1, 2}), f.graph->qprevs[3]);
}

TEST(GraphCtor, TruthVariantPairedWithBypass) {
    // Each truth variant node is immediately followed by a reference-allele bypass node over the
    // same reference span, carrying tskips == the bypassed variant index. calc_prec_recall relies
    // on that adjacency: excise_bypass() reads the variant node at bypass_tni - 1.
    GraphFixture f = build_fixture("ACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    ASSERT_EQ(4, f.graph->tnodes);
    EXPECT_EQ(std::vector<std::string>({"_C", "_A", "_G", "_TA"}), f.graph->tseqs);
    EXPECT_EQ(std::vector<int>({0, 1, 1, 2}), f.graph->tbegs);
    EXPECT_EQ(std::vector<int>({1, 2, 2, 4}), f.graph->tends);
    EXPECT_EQ(std::vector<int>({TYPE_REF, TYPE_SUB, TYPE_REF, TYPE_REF}), f.graph->ttypes);
    EXPECT_EQ(std::vector<int>({-1, 0, -1, -1}), f.graph->tidxs);
    EXPECT_EQ(std::vector<int>({-1, -1, 0, -1}), f.graph->tskips);
}

TEST(GraphCtor, TruthStringExcludesBypassAllele) {
    // this->truth is the selected truth haplotype: the alt is spliced in and the bypass node's
    // reference allele is absent, even though both nodes exist in the graph.
    GraphFixture f = build_fixture("ACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    EXPECT_EQ("CATA", f.graph->truth);
    EXPECT_EQ("CGTA", f.graph->ref) << "the reference keeps the bypassed allele";
}

TEST(GraphCtor, TruthHapFilterExcludesOtherHap) {
    // GT_REF_ALT1 places the truth variant on haplotype 1 only, so a HAP1 (index 0) graph skips
    // it entirely and collapses to one reference node, while a HAP2 graph builds alt and bypass.
    const std::vector<var_desc> qv = {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}};
    const std::vector<var_desc> tv = {{2, 1, TYPE_SUB, "G", "A", GT_REF_ALT1, 60, 0, 0}};

    GraphFixture hap1 = build_fixture("ACGTACGT", qv, tv, HAP1);
    EXPECT_EQ(1, hap1.graph->tnodes);
    EXPECT_EQ("CGTA", hap1.graph->truth) << "hap1 truth should be pure reference";

    GraphFixture hap2 = build_fixture("ACGTACGT", qv, tv, HAP2);
    EXPECT_EQ(4, hap2.graph->tnodes);
    EXPECT_EQ("CATA", hap2.graph->truth);
}

TEST(GraphCtor, SkipEvaluatedQueryVariantOmitted) {
    // A query variant already labeled on this haplotype (errtype != ERRTYPE_UN) contributes no
    // node, so only reference nodes remain. Its position still sets the graph's span.
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto qvars = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 60, 0, 0}});
    qvars->errtypes[HAP1][0] = ERRTYPE_TP;
    auto sc = make_ctgSuperclusters(qvars, make_ctgVariants("chr1", {}));

    auto graph = make_graph(sc, ref, "chr1", HAP1);
    ASSERT_EQ(1, graph->qnodes);
    EXPECT_EQ(TYPE_REF, graph->qtypes[0]);
    EXPECT_EQ("_CGTA", graph->qseqs[0]);
}

TEST(GraphCtor, SkipEvaluatedTruthVariantOmitted) {
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto tvars = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});
    tvars->errtypes[HAP1][0] = ERRTYPE_FN;
    auto sc = make_ctgSuperclusters(make_ctgVariants("chr1", {}), tvars);

    auto graph = make_graph(sc, ref, "chr1", HAP1);
    ASSERT_EQ(1, graph->tnodes);
    EXPECT_EQ(-1, graph->tskips[0]) << "no bypass node for a skipped truth variant";
    EXPECT_EQ("CGTA", graph->truth);
}

TEST(GraphCtor, TruthLinearChain) {
    // Two truth SNPs at pos 2 and 6 over "ACGTACGTACGT": the truth side is a linear chain of
    // (ref, alt, bypass) triples closed by a trailing reference node, spanning [1, 8].
    GraphFixture f = build_fixture("ACGTACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
             {6, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    ASSERT_EQ(7, f.graph->tnodes);
    EXPECT_EQ(std::vector<std::string>({"_C", "_A", "_G", "_TAC", "_A", "_G", "_TA"}),
            f.graph->tseqs);
    // the bypass nodes are 2 and 5, each naming the variant it routes around
    EXPECT_EQ(std::vector<int>({-1, -1, 0, -1, -1, 1, -1}), f.graph->tskips);
    EXPECT_EQ(std::vector<int>({-1, 0, -1, -1, 1, -1, -1}), f.graph->tidxs);
    EXPECT_EQ("CATACATA", f.graph->truth);
    EXPECT_EQ("CGTACGTA", f.graph->ref);

    // each triple's alt and bypass are parallel: both lead to the following reference node
    EXPECT_EQ(std::vector<int>({1, 2}), f.graph->tnexts[0]);
    EXPECT_EQ(std::vector<int>({3}), f.graph->tnexts[1]);
    EXPECT_EQ(std::vector<int>({3}), f.graph->tnexts[2]);
    EXPECT_EQ(std::vector<int>({4, 5}), f.graph->tnexts[3]);
    EXPECT_EQ(std::vector<int>({6}), f.graph->tnexts[4]);
    EXPECT_EQ(std::vector<int>({6}), f.graph->tnexts[5]);
}

TEST(GraphCtor, InsertionNodeIsZeroWidth) {
    // An insertion consumes no reference, so its node has qbegs == qends and pushes nothing onto
    // qnode_ends -- there is no parallel reference-allele node, only a direct ref-to-ref edge.
    GraphFixture f = build_fixture("ACGTACGT",
            {{3, 0, TYPE_INS, "", "TT", GT_ALT1_REF, 60, 0, 0}}, {});

    ASSERT_EQ(3, f.graph->qnodes);
    EXPECT_EQ("_TT", f.graph->qseqs[1]);
    EXPECT_EQ(TYPE_INS, f.graph->qtypes[1]);
    EXPECT_EQ(f.graph->qbegs[1], f.graph->qends[1]);
    // node 0 leads into both the insertion and the following reference node
    EXPECT_EQ(std::vector<int>({1, 2}), f.graph->qnexts[0]);
}

TEST(GraphCtor, DeletionNodeSpansRefAllele) {
    // A deletion emits an empty alt but still spans its two reference bases, and its end is
    // pushed onto qnode_ends, so the deleted reference reappears as a parallel node.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 2, TYPE_DEL, "GT", "", GT_ALT1_REF, 60, 0, 0}}, {});

    ASSERT_EQ(4, f.graph->qnodes);
    EXPECT_EQ("_", f.graph->qseqs[1]);
    EXPECT_EQ(TYPE_DEL, f.graph->qtypes[1]);
    EXPECT_EQ(2, f.graph->qends[1] - f.graph->qbegs[1]);
    EXPECT_EQ("_GT", f.graph->qseqs[2]) << "the deleted reference allele runs in parallel";
    EXPECT_EQ(f.graph->qbegs[1], f.graph->qbegs[2]);
    EXPECT_EQ(f.graph->qends[1], f.graph->qends[2]);
}

TEST(GraphCtor, OverlappingQueryVariantsShareRefSpan) {
    // A 3bp deletion at pos 2 spans the SNP at pos 3. qnode_ends orders the cuts, so the graph
    // holds both a deletion node leaping [1,4) and a SNP path through it.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 3, TYPE_DEL, "GTA", "", GT_ALT1_REF, 60, 0, 0},
             {3, 1, TYPE_SUB, "T", "C", GT_ALT1_REF, 60, 0, 0}}, {});

    ASSERT_EQ(7, f.graph->qnodes);
    EXPECT_EQ(std::vector<std::string>({"_C", "_", "_G", "_C", "_T", "_A", "_CG"}), f.graph->qseqs);
    EXPECT_EQ(std::vector<int>({0, 1, 1, 2, 2, 3, 4}), f.graph->qbegs);
    EXPECT_EQ(std::vector<int>({1, 4, 2, 3, 3, 4, 6}), f.graph->qends);
    EXPECT_EQ(std::vector<int>({-1, 0, -1, 1, -1, -1, -1}), f.graph->qidxs);

    // the deletion node jumps straight from offset 1 to offset 4, skipping the SNP
    EXPECT_EQ(std::vector<int>({6}), f.graph->qnexts[1]);
    // while the SNP is reached through the deletion's parallel reference allele
    EXPECT_EQ(std::vector<int>({3, 4}), f.graph->qnexts[2]);
}

// A truth insertion is a zero-width (in reference coordinates) locus. When another truth
// variant abuts it, a direct edge that leaps the insertion must not exist: every path across
// the locus has to route through the insertion's own alt or bypass node. This asserts that the
// reference node immediately after the insertion has only zero-width predecessors, so no
// reference-spanning node (e.g. the neighbouring SUB's alt or bypass) can bypass the insertion
// for free. Without this the insertion would be silently skipped and left unlabeled.
TEST(GraphInsertionEdges, AdjacentVariantCannotLeapInsertion) {
    auto ref = make_fasta("chr1", "ACGTACGTAC");

    // truth SUB at pos 2 (G->T) immediately followed by truth INS at pos 3 (->TTT), both hap0
    auto tv = make_ctgVariants("chr1", {
            {2, 1, TYPE_SUB, "G", "T",   GT_ALT1_REF, 60, 0, 0},
            {3, 0, TYPE_INS, "",  "TTT", GT_ALT1_REF, 60, 0, 0}});
    auto sc = make_ctgSuperclusters(make_ctgVariants("chr1", {}), tv);

    auto graph = make_graph(sc, ref, "chr1", HAP1);

    // locate the insertion locus (the zero-width truth variant node)
    int ins_coord = -1;
    for (int tn = 0; tn < graph->tnodes; tn++)
        if (graph->tidxs[tn] >= 0 && graph->tbegs[tn] == graph->tends[tn])
            ins_coord = graph->tbegs[tn];
    ASSERT_GE(ins_coord, 0) << "no zero-width insertion node found";

    // locate the reference node immediately to the right of the insertion locus
    int right_ref = -1;
    for (int tn = 0; tn < graph->tnodes; tn++)
        if (graph->ttypes[tn] == TYPE_REF && graph->tskips[tn] < 0 && graph->tidxs[tn] < 0 &&
                graph->tbegs[tn] == ins_coord && graph->tends[tn] > graph->tbegs[tn])
            right_ref = tn;
    ASSERT_GE(right_ref, 0) << "no reference node found after the insertion locus";

    // every predecessor must be a zero-width node at the locus (the insertion's alt/bypass);
    // a reference-spanning predecessor would be a free leap over the insertion
    ASSERT_FALSE(graph->tprevs[right_ref].empty());
    for (int p : graph->tprevs[right_ref])
        EXPECT_EQ(graph->tbegs[p], graph->tends[p])
            << "node " << p << " (spanning " << graph->tbegs[p] << ".." << graph->tends[p]
            << ") leaps the insertion at " << ins_coord;
}

/* Graph::get_truth_pos ***************************************************************************/

// get_truth_pos returns the number of this->truth characters consumed on arriving at cell
// (tn, ti), so it is the exclusive upper bound of the consumed prefix: the character the cell
// itself matches is truth[get_truth_pos(tn, ti) - 1]. Each tseq starts with a '_' placeholder,
// which is why every node contributes size() - 1.

TEST(GetTruthPos, FirstNodeIsIdentity) {
    // nothing precedes node 0, so the offset within the node is the answer
    GraphFixture f = build_fixture("ACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    EXPECT_EQ(0, f.graph->get_truth_pos(0, 0));
    EXPECT_EQ(1, f.graph->get_truth_pos(0, 1));
}

TEST(GetTruthPos, SkipsBypassNodes) {
    // nodes are (ref "_C", alt "_A", bypass "_G", ref "_TA"). Reaching node 3 consumes "CA", so
    // the answer is 2. Counting the bypass node's single base as well would give 3.
    GraphFixture f = build_fixture("ACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    ASSERT_EQ(4, f.graph->tnodes);
    ASSERT_EQ(0, f.graph->tskips[2]);
    EXPECT_EQ(2, f.graph->get_truth_pos(3, 0));
    EXPECT_EQ(4, f.graph->get_truth_pos(3, 2));
    EXPECT_EQ(int(f.graph->truth.size()), f.graph->get_truth_pos(3, 2));
}

TEST(GetTruthPos, SubtractsUnderscorePerNode) {
    // seven nodes "_C", "_A", bypass "_G", "_TAC", "_A", bypass "_G", "_TA" over truth
    // "CATACATA": reaching node 6 consumes 1 + 1 + 3 + 1 = 6 characters, the two bypasses
    // contributing nothing.
    GraphFixture f = build_fixture("ACGTACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
             {6, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    ASSERT_EQ(7, f.graph->tnodes);
    EXPECT_EQ(1, f.graph->get_truth_pos(1, 0));
    EXPECT_EQ(2, f.graph->get_truth_pos(3, 0));
    EXPECT_EQ(5, f.graph->get_truth_pos(4, 0));
    EXPECT_EQ(6, f.graph->get_truth_pos(6, 0));
    EXPECT_EQ(8, f.graph->get_truth_pos(6, 2));
}

TEST(GetTruthPos, CountsTruthCharactersConsumed) {
    // The invariant every caller depends on: for every non-bypass node and every offset past the
    // '_' placeholder, the cell's own character is the last one the returned prefix contains.
    GraphFixture f = build_fixture("ACGTACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
             {6, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    for (int tn = 0; tn < f.graph->tnodes; tn++) {
        if (f.graph->tskips[tn] >= 0) continue; // not part of this->truth
        for (int ti = 1; ti < int(f.graph->tseqs[tn].size()); ti++) {
            const int pos = f.graph->get_truth_pos(tn, ti);
            ASSERT_GE(pos, 1) << "node " << tn << " offset " << ti;
            ASSERT_LE(pos, int(f.graph->truth.size())) << "node " << tn << " offset " << ti;
            EXPECT_EQ(f.graph->tseqs[tn][ti], f.graph->truth[pos - 1])
                << "node " << tn << " offset " << ti;
        }
    }
}

/* calc_prec_recall_aln ***************************************************************************/

// Cost is unit per edit, plus a fractional (1 - g.credit_threshold) * len toll for entering a
// truth variant's reference-allele bypass node, where len is the longer of the variant's ref and
// alt alleles. The return value is the total rounded up: int(ceil(cost - EPSILON)).

TEST(PrecRecallAln, IdenticalToReferenceScoresZero) {
    GlobalsGuard guard;

    // the query SNP's parallel reference allele reproduces the truth exactly, at no cost
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 60, 0, 0}}, {});

    std::unordered_map<idx4, idx4> ptrs;
    EXPECT_EQ(0, calc_prec_recall_aln(f.graph, ptrs, false));
}

TEST(PrecRecallAln, MatchingVariantCostsZero) {
    GlobalsGuard guard;

    // query "CATA" through its alt equals truth "CATA" through its alt: no edits, no toll
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    ASSERT_EQ("CATA", f.graph->truth);
    std::unordered_map<idx4, idx4> ptrs;
    EXPECT_EQ(0, calc_prec_recall_aln(f.graph, ptrs, false));
}

TEST(PrecRecallAln, MissedSnpCostsFractionalBypassToll) {
    GlobalsGuard guard;

    // The query is pure reference, so the only route past the truth SNP is its bypass node:
    // (1 - 0.98) * max(|"G"|, |"A"|) = 0.02, which ceil()s to 1. Reproducing the SNP instead
    // would cost a full substitution.
    GraphFixture f = build_fixture("ACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    ASSERT_DOUBLE_EQ(0.98, g.credit_threshold);
    std::unordered_map<idx4, idx4> ptrs;
    EXPECT_EQ(1, calc_prec_recall_aln(f.graph, ptrs, false));
}

TEST(PrecRecallAln, BypassTollScalesWithVariantLength) {
    GlobalsGuard guard;

    // A missed 60bp deletion tolls (1 - 0.98) * 60 = 1.2, which ceil()s to 2. A length-independent
    // toll, or truncation instead of ceil(), would report 1.
    const std::string ref_seq = periodic_ref(70);
    GraphFixture f = build_fixture(ref_seq, {},
            {{5, 60, TYPE_DEL, ref_seq.substr(5, 60), "", GT_ALT1_REF, 60, 0, 0}});

    std::unordered_map<idx4, idx4> ptrs;
    EXPECT_EQ(2, calc_prec_recall_aln(f.graph, ptrs, false));
}

TEST(PrecRecallAln, BypassTollsAccumulate) {
    GlobalsGuard guard;

    // two missed 1bp SNPs toll 0.02 each; 0.04 still ceil()s to 1
    GraphFixture f = build_fixture("ACGTACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
             {6, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    std::unordered_map<idx4, idx4> ptrs;
    EXPECT_EQ(1, calc_prec_recall_aln(f.graph, ptrs, false));
}

TEST(PrecRecallAln, EditsCostOnePerBase) {
    GlobalsGuard guard;

    // With credit_threshold 0 the bypass toll is the full 1.0 * len = 3, so the cheaper route is
    // through the alt alleles and the score is the plain edit count between them. Truth alt "GGG"
    // against query alt "GG" is one deleted base, and against "G" it is two -- there is no
    // gap-open discount, unlike wf_swg_align's affine scoring.
    g.credit_threshold = 0;

    GraphFixture one = build_fixture("ACGTACGTACGT",
            {{3, 3, TYPE_CPX, "TAC", "GG", GT_ALT1_REF, 60, 0, 0}},
            {{3, 3, TYPE_CPX, "TAC", "GGG", GT_ALT1_REF, 60, 0, 0}});
    std::unordered_map<idx4, idx4> one_ptrs;
    EXPECT_EQ(1, calc_prec_recall_aln(one.graph, one_ptrs, false));

    GraphFixture two = build_fixture("ACGTACGTACGT",
            {{3, 3, TYPE_CPX, "TAC", "G", GT_ALT1_REF, 60, 0, 0}},
            {{3, 3, TYPE_CPX, "TAC", "GGG", GT_ALT1_REF, 60, 0, 0}});
    std::unordered_map<idx4, idx4> two_ptrs;
    EXPECT_EQ(2, calc_prec_recall_aln(two.graph, two_ptrs, false));
}

TEST(PrecRecallAln, PtrsTraceEndpointBackToStart) {
    GlobalsGuard guard;

    // calc_prec_recall walks ptrs from the endpoint until it reaches the idx4(0,0,-1,-1) sentinel,
    // so every cell on that chain must be present and the walk must terminate.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    std::unordered_map<idx4, idx4> ptrs;
    calc_prec_recall_aln(f.graph, ptrs, false);

    const idx4 end(f.graph->qnodes-1, f.graph->tnodes-1,
            int(f.graph->qseqs[f.graph->qnodes-1].size()) - 1,
            int(f.graph->tseqs[f.graph->tnodes-1].size()) - 1);
    ASSERT_EQ(size_t(1), ptrs.count(end)) << "endpoint missing from the pointer map";

    const idx4 sentinel(0, 0, -1, -1);
    idx4 curr = end;
    int steps = 0;
    while (curr != sentinel) {
        ASSERT_EQ(size_t(1), ptrs.count(curr)) << "chain broke after " << steps << " steps";
        curr = ptrs.at(curr);
        ASSERT_LT(++steps, 1000) << "pointer chain does not terminate";
    }
    EXPECT_GT(steps, 1) << "the chain should cross more than the start cell";
}

TEST(PrecRecallAln, UnreachableEndpointErrors) {
    GlobalsGuard guard;

    // No make_graph() input leaves the endpoint unreachable, so the successor lists are cleared by
    // hand: with no way out of node 0 the queue drains before the endpoint and the guard fires.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 60, 0, 0}}, {});
    ASSERT_GT(f.graph->qnodes, 1) << "the endpoint must lie in a later node to be unreachable";
    for (std::vector<int> & nexts : f.graph->qnexts) nexts.clear();
    for (std::vector<int> & nexts : f.graph->tnexts) nexts.clear();

    std::unordered_map<idx4, idx4> ptrs;
    EXPECT_EXIT(calc_prec_recall_aln(f.graph, ptrs, false),
            testing::ExitedWithCode(1), "Endpoint unreachable");
}

/* calc_prec_recall *******************************************************************************/

TEST(PrecRecall, UnmatchedQueryVariantDefaultsToFalsePositive) {
    GlobalsGuard guard;

    // With no truth variants the query SNP's parallel reference allele wins the alignment, so the
    // SNP never joins a sync group and keeps the default FP the labeling pass seeds it with. Its
    // callq is its own variant quality, not a group minimum.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 30, 0, 0}}, {});

    align_and_label(f.graph, HAP1);

    EXPECT_EQ(ERRTYPE_FP, f.qvars->errtypes[HAP1][0]);
    EXPECT_FLOAT_EQ(30.0f, f.qvars->callq[HAP1][0]);
}

TEST(PrecRecall, ReproducedVariantIsTruePositive) {
    GlobalsGuard guard;

    // Query and truth both carry the pos-2 SNP. The group spans reference "GTA" against truth
    // "ATA", one edit apart, and the query needed no edits: credit = 1 - 0/1 = 1 >= 0.98.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    align_and_label(f.graph, HAP1);

    EXPECT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_TP, f.tvars->errtypes[HAP1][0]);
    EXPECT_FLOAT_EQ(1.0f, f.qvars->credit[HAP1][0]);
    EXPECT_EQ(1, f.qvars->ref_ed[HAP1][0]);
    EXPECT_EQ(0, f.qvars->query_ed[HAP1][0]);
    EXPECT_EQ(0, f.qvars->sync_group[HAP1][0]);
    EXPECT_FLOAT_EQ(60.0f, f.qvars->callq[HAP1][0]);
    // a TP sets the query variant's calculated genotype on the matched haplotype
    EXPECT_TRUE(f.qvars->var_on_hap(0, HAP1, true /* calc */));
}

TEST(PrecRecall, MissedVariantIsFalseNegativeViaBypass) {
    GlobalsGuard guard;

    // A pure-reference query routes around the truth SNP through its bypass node. The FN is
    // labeled at that transition, outside any sync group, so its credit and edit distances are
    // zeroed rather than measured.
    GraphFixture f = build_fixture("ACGTACGT", {},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    align_and_label(f.graph, HAP1);

    EXPECT_EQ(ERRTYPE_FN, f.tvars->errtypes[HAP1][0]);
    EXPECT_FLOAT_EQ(0.0f, f.tvars->credit[HAP1][0]);
    EXPECT_EQ(0, f.tvars->ref_ed[HAP1][0]);
    EXPECT_EQ(0, f.tvars->query_ed[HAP1][0]);
    EXPECT_EQ(0, f.tvars->sync_group[HAP1][0]);
}

// End-to-end guard: with the leap-suppression edge rule, a truth insertion abutting another
// truth variant is forced onto the alignment path (via its alt or bypass node) and is labeled by
// the backtrack itself -- not by the removed safety sweep. Here the query equals the reference, so
// both truth variants are missed and must be labeled FN via their bypass nodes, never left UNKNOWN.
TEST(GraphInsertionEdges, AdjacentInsertionLabeledWithoutSweep) {
    auto ref = make_fasta("chr1", "ACGTACGTAC");

    auto tv = make_ctgVariants("chr1", {
            {2, 1, TYPE_SUB, "G", "T",   GT_ALT1_REF, 60, 0, 0},
            {3, 0, TYPE_INS, "",  "TTT", GT_ALT1_REF, 60, 0, 0}});
    auto sc = make_ctgSuperclusters(make_ctgVariants("chr1", {}), tv);

    auto graph = make_graph(sc, ref, "chr1", HAP1);
    std::unordered_map<idx4, idx4> ptrs;
    calc_prec_recall_aln(graph, ptrs, false);
    calc_prec_recall(graph, ptrs, HAP1, false);

    // query == reference: neither truth variant is reproduced, both are FN (not left UNKNOWN)
    EXPECT_EQ(ERRTYPE_FN, tv->errtypes[HAP1][0]) << "SUB should be a false negative";
    EXPECT_EQ(ERRTYPE_FN, tv->errtypes[HAP1][1]) << "insertion should be a false negative";
}

// Regression guard for consecutive bypassed (FN) truth variants. Two adjacent truth SNPs the query
// misses are both routed through their reference-allele bypass nodes (FN). Each bypass span must be
// excised from the enclosing sync group's credit measurement -- not just the first one reached
// during backtracking. A neighbouring reproduced SNP (TP) shares the same group; with only one span
// excised its reference edit distance is inflated (2 instead of 1), so ref_ed pins the excision of
// BOTH bypasses.
TEST(GraphBypass, ConsecutiveBypassesBothExcised) {
    auto ref = make_fasta("chr1", "ACGTACGT");

    // truth: SNP pos2 (reproduced -> TP), then adjacent SNPs pos3 and pos4 (both missed -> FN)
    auto tv = make_ctgVariants("chr1", {
            {2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
            {3, 1, TYPE_SUB, "T", "G", GT_ALT1_REF, 60, 0, 0},
            {4, 1, TYPE_SUB, "A", "C", GT_ALT1_REF, 60, 0, 0}});

    // query: reproduces only the pos2 SNP, so pos3 and pos4 are missed
    auto qv = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    auto sc = make_ctgSuperclusters(qv, tv);
    auto graph = make_graph(sc, ref, "chr1", HAP1);
    std::unordered_map<idx4, idx4> ptrs;
    calc_prec_recall_aln(graph, ptrs, false);
    calc_prec_recall(graph, ptrs, HAP1, false);

    EXPECT_EQ(ERRTYPE_TP, tv->errtypes[HAP1][0]) << "reproduced pos2 SNP should be TP";
    EXPECT_EQ(ERRTYPE_FN, tv->errtypes[HAP1][1]) << "missed pos3 SNP should be FN";
    EXPECT_EQ(ERRTYPE_FN, tv->errtypes[HAP1][2]) << "missed pos4 SNP should be FN";
    // both bypass spans excised: the TP group's ref edit distance counts only the pos2 SNP.
    // Leaving the second bypass un-excised inflates this to 2.
    EXPECT_EQ(1, tv->ref_ed[HAP1][0]) << "consecutive bypass span not fully excised";
}

TEST(PrecRecall, ImperfectMatchLosesToBypassAtDefaultThreshold) {
    GlobalsGuard guard;

    // A consequence of the default credit_threshold that is easy to overlook: the bypass toll for
    // a 2bp truth variant is only (1 - 0.98) * 2 = 0.04, far below the cost of a single edit, so
    // a query that reproduces the variant imperfectly loses to the bypass. The truth variant is
    // labeled FN rather than earning partial credit, and the query call keeps its default FP.
    GraphFixture f = build_fixture("ACGTACGTACGT",
            {{3, 1, TYPE_SUB, "T", "G", GT_ALT1_REF, 60, 0, 0}},
            {{3, 2, TYPE_CPX, "TA", "GG", GT_ALT1_REF, 60, 0, 0}});

    align_and_label(f.graph, HAP1);

    EXPECT_EQ(ERRTYPE_FN, f.tvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_FP, f.qvars->errtypes[HAP1][0]);
    EXPECT_FLOAT_EQ(0.0f, f.tvars->credit[HAP1][0]);
}

/**
 * @brief Aligns the half-reproduced 2bp truth variant with a toll high enough to force the alt.
 *
 * At credit_threshold 0.4 the bypass tolls 0.6 * 2 = 1.2, above the single edit the alt route
 * costs, so the alignment runs through the truth alt and the group measures reference "TACG"
 * against truth "GGCG" (two edits) with one query edit: credit = 1 - 1/2 = 0.5. Labeling is left
 * to the caller, which sets its own threshold to test the comparison in emit_sync_group.
 */
static GraphFixture half_credit_alignment(std::unordered_map<idx4, idx4> & ptrs) {
    g.credit_threshold = 0.4;
    GraphFixture f = build_fixture("ACGTACGTACGT",
            {{3, 1, TYPE_SUB, "T", "G", GT_ALT1_REF, 60, 0, 0}},
            {{3, 2, TYPE_CPX, "TA", "GG", GT_ALT1_REF, 60, 0, 0}});
    calc_prec_recall_aln(f.graph, ptrs, false);
    return f;
}

TEST(PrecRecall, CreditAtThresholdIsTruePositive) {
    GlobalsGuard guard;

    std::unordered_map<idx4, idx4> ptrs;
    GraphFixture f = half_credit_alignment(ptrs);

    // the comparison is >=, so credit exactly at the threshold is a true positive
    g.credit_threshold = 0.5;
    calc_prec_recall(f.graph, ptrs, HAP1, false);

    EXPECT_FLOAT_EQ(0.5f, f.qvars->credit[HAP1][0]);
    EXPECT_EQ(2, f.qvars->ref_ed[HAP1][0]);
    EXPECT_EQ(1, f.qvars->query_ed[HAP1][0]);
    EXPECT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_TP, f.tvars->errtypes[HAP1][0]);
}

TEST(PrecRecall, CreditBelowThresholdIsFalsePositiveAndFalseNegative) {
    GlobalsGuard guard;

    std::unordered_map<idx4, idx4> ptrs;
    GraphFixture f = half_credit_alignment(ptrs);

    // one point above the same credit, and the group fails: the query call is a false positive
    // and the truth variant in the same group is downgraded to a false negative
    g.credit_threshold = 0.51;
    calc_prec_recall(f.graph, ptrs, HAP1, false);

    EXPECT_FLOAT_EQ(0.5f, f.qvars->credit[HAP1][0]);
    EXPECT_EQ(ERRTYPE_FP, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_FN, f.tvars->errtypes[HAP1][0]);
}

TEST(PrecRecall, ZeroReferenceDistanceGivesZeroCredit) {
    GlobalsGuard guard;

    // emit_sync_group() guards against dividing by a zero reference distance. No make_graph()
    // fixture reaches that guard -- a query alt only wins the alignment where the truth differs
    // from the reference, which forces ref_dist >= 1 -- so graph->ref is overwritten with the
    // equal-length graph->truth after aligning, leaving the two spliced spans identical.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    std::unordered_map<idx4, idx4> ptrs;
    calc_prec_recall_aln(f.graph, ptrs, false);
    ASSERT_EQ(f.graph->ref.size(), f.graph->truth.size());
    f.graph->ref = f.graph->truth;
    calc_prec_recall(f.graph, ptrs, HAP1, false);

    EXPECT_EQ(0, f.qvars->ref_ed[HAP1][0]);
    EXPECT_FLOAT_EQ(0.0f, f.qvars->credit[HAP1][0]);
    // credit 0 is below the threshold, so the otherwise-perfect call becomes FP and FN
    EXPECT_EQ(ERRTYPE_FP, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_FN, f.tvars->errtypes[HAP1][0]);
}

TEST(PrecRecall, CallQualityIsGroupMinimum) {
    GlobalsGuard guard;

    // Two query SNPs at pos 3 and 4 reproduce one 2bp truth variant. Moving between the two query
    // alt nodes stays inside the single truth node, so it is neither a same-submatrix reference
    // move nor a change of both submatrices: no sync point separates them and both land in one
    // group, whose call quality is the minimum over its members (min(30, 50) = 30).
    GraphFixture f = build_fixture("ACGTACGTACGT",
            {{3, 1, TYPE_SUB, "T", "G", GT_ALT1_REF, 30, 0, 0},
             {4, 1, TYPE_SUB, "A", "G", GT_ALT1_REF, 50, 0, 0}},
            {{3, 2, TYPE_CPX, "TA", "GG", GT_ALT1_REF, 60, 0, 0}});

    align_and_label(f.graph, HAP1);

    ASSERT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][0]);
    ASSERT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][1]);
    EXPECT_EQ(f.qvars->sync_group[HAP1][0], f.qvars->sync_group[HAP1][1])
        << "both calls should share one sync group";
    EXPECT_FLOAT_EQ(30.0f, f.qvars->callq[HAP1][0]);
    EXPECT_FLOAT_EQ(30.0f, f.qvars->callq[HAP1][1]);
    EXPECT_FLOAT_EQ(30.0f, f.tvars->callq[HAP1][0]) << "the truth variant shares the group quality";
    // reference "TACG" against truth "GGCG" is two edits, both of which the query reproduced
    EXPECT_EQ(2, f.qvars->ref_ed[HAP1][0]);
    EXPECT_EQ(0, f.qvars->query_ed[HAP1][0]);
}

TEST(PrecRecall, SyncGroupIncrementsPerGroup) {
    GlobalsGuard guard;

    // Two reproduced SNPs far enough apart to sync between them get their own groups. The
    // backtrack runs right to left, so the rightmost group is numbered 0.
    GraphFixture f = build_fixture("ACGTACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
             {6, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
             {6, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    align_and_label(f.graph, HAP1);

    EXPECT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][1]);
    EXPECT_EQ(1, f.qvars->sync_group[HAP1][0]);
    EXPECT_EQ(0, f.qvars->sync_group[HAP1][1]);
    EXPECT_EQ(1, f.tvars->sync_group[HAP1][0]);
    EXPECT_EQ(0, f.tvars->sync_group[HAP1][1]);
    // each group measures only its own SNP, so one edit apart from the reference
    EXPECT_EQ(1, f.qvars->ref_ed[HAP1][0]);
    EXPECT_EQ(1, f.qvars->ref_ed[HAP1][1]);
}

/* evaluate_variants ******************************************************************************/

TEST(EvaluateVariants, MatchesManualAlignAndLabel) {
    GlobalsGuard guard;

    // evaluate_variants is exactly graph construction followed by one alignment and one labeling
    // pass, so it must agree with that sequence run by hand on the same input.
    const std::vector<var_desc> qv = {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
                                      {6, 1, TYPE_SUB, "G", "T", GT_ALT1_REF, 40, 0, 0}};
    const std::vector<var_desc> tv = {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0},
                                      {6, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}};

    GraphFixture manual = build_fixture("ACGTACGTACGT", qv, tv);
    align_and_label(manual.graph, HAP1);

    GraphFixture via_wrapper = build_fixture("ACGTACGTACGT", qv, tv);
    evaluate_variants(via_wrapper.sc, 0, via_wrapper.ref, "chr1", HAP1, false);

    for (int i = 0; i < 2; i++) {
        EXPECT_EQ(manual.qvars->errtypes[HAP1][i], via_wrapper.qvars->errtypes[HAP1][i]) << i;
        EXPECT_EQ(manual.qvars->sync_group[HAP1][i], via_wrapper.qvars->sync_group[HAP1][i]) << i;
        EXPECT_EQ(manual.qvars->ref_ed[HAP1][i], via_wrapper.qvars->ref_ed[HAP1][i]) << i;
        EXPECT_EQ(manual.qvars->query_ed[HAP1][i], via_wrapper.qvars->query_ed[HAP1][i]) << i;
        EXPECT_FLOAT_EQ(manual.qvars->credit[HAP1][i], via_wrapper.qvars->credit[HAP1][i]) << i;
        EXPECT_FLOAT_EQ(manual.qvars->callq[HAP1][i], via_wrapper.qvars->callq[HAP1][i]) << i;
        EXPECT_EQ(manual.tvars->errtypes[HAP1][i], via_wrapper.tvars->errtypes[HAP1][i]) << i;
        EXPECT_FLOAT_EQ(manual.tvars->credit[HAP1][i], via_wrapper.tvars->credit[HAP1][i]) << i;
    }
}

TEST(EvaluateVariants, LabelsEachHaplotypeIndependently) {
    GlobalsGuard guard;

    // The truth variant is on haplotype 0 only. Evaluating haplotype 0 pairs it with the query
    // call; evaluating haplotype 1 sees a pure-reference truth, so it never labels the truth
    // variant at all and the query call stays FP. Neither pass touches the other's lane.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});

    evaluate_variants(f.sc, 0, f.ref, "chr1", HAP1, false);
    EXPECT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_TP, f.tvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_UN, f.qvars->errtypes[HAP2][0]) << "hap2 lane written by the hap1 pass";
    EXPECT_EQ(ERRTYPE_UN, f.tvars->errtypes[HAP2][0]) << "hap2 lane written by the hap1 pass";

    evaluate_variants(f.sc, 0, f.ref, "chr1", HAP2, false);
    EXPECT_EQ(ERRTYPE_FP, f.qvars->errtypes[HAP2][0]);
    EXPECT_EQ(ERRTYPE_UN, f.tvars->errtypes[HAP2][0]) << "off-haplotype truth is never evaluated";
    EXPECT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][0]) << "hap1 lane clobbered by the hap2 pass";
}

TEST(EvaluateVariants, PreservesAlreadyEvaluatedVariants) {
    GlobalsGuard guard;

    // A query variant labeled before the call is excluded from the graph, so neither the
    // default-FP seeding nor the backtrack can overwrite its label or its call quality.
    auto ref = make_fasta("chr1", "ACGTACGT");
    auto qvars = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});
    auto tvars = make_ctgVariants("chr1", {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});
    qvars->superclusters[0] = 0;
    tvars->superclusters[0] = 0;
    qvars->errtypes[HAP1][0] = ERRTYPE_TP;
    qvars->callq[HAP1][0] = 17;
    auto sc = make_ctgSuperclusters(qvars, tvars);

    evaluate_variants(sc, 0, ref, "chr1", HAP1, false);

    EXPECT_EQ(ERRTYPE_TP, qvars->errtypes[HAP1][0]);
    EXPECT_FLOAT_EQ(17.0f, qvars->callq[HAP1][0]);
    // the truth variant is still evaluated, and with no query allele left it is a false negative
    EXPECT_EQ(ERRTYPE_FN, tvars->errtypes[HAP1][0]);
}

TEST(EvaluateVariants, PrintDoesNotChangeLabels) {
    GlobalsGuard guard;

    // print=true adds Graph::print() and the backtrack trace to stdout without affecting any
    // label. Stdout is captured so the trace does not bury the test log.
    const std::vector<var_desc> qv = {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}};
    const std::vector<var_desc> tv = {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}};

    GraphFixture quiet = build_fixture("ACGTACGT", qv, tv);
    evaluate_variants(quiet.sc, 0, quiet.ref, "chr1", HAP1, false);

    GraphFixture loud = build_fixture("ACGTACGT", qv, tv);
    testing::internal::CaptureStdout();
    evaluate_variants(loud.sc, 0, loud.ref, "chr1", HAP1, true);
    const std::string output = testing::internal::GetCapturedStdout();

    EXPECT_FALSE(output.empty()) << "print=true should emit a trace";
    EXPECT_EQ(quiet.qvars->errtypes[HAP1][0], loud.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(quiet.tvars->errtypes[HAP1][0], loud.tvars->errtypes[HAP1][0]);
    EXPECT_EQ(quiet.qvars->ref_ed[HAP1][0], loud.qvars->ref_ed[HAP1][0]);
    EXPECT_FLOAT_EQ(quiet.qvars->credit[HAP1][0], loud.qvars->credit[HAP1][0]);
}

/* precision_recall_wrapper ***********************************************************************/

TEST(PrecisionRecallWrapper, EmptyRangeReturns) {
    GlobalsGuard guard;

    // stop == start returns before sc_groups is ever indexed, so an empty grouping is safe and
    // no variant is touched.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});
    auto sc_data = make_superclusterData({"chr1"}, {8}, {2}, {f.sc}, f.ref);

    const std::vector< std::vector< std::vector<int> > > sc_groups;
    precision_recall_wrapper(sc_data.get(), sc_groups, 0, 0, 0, false, false);

    EXPECT_EQ(ERRTYPE_UN, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_UN, f.tvars->errtypes[HAP1][0]);
}

TEST(PrecisionRecallWrapper, EvaluatesBothHaplotypesForOneSupercluster) {
    GlobalsGuard guard;

    // One supercluster on contig 0 is evaluated for both haplotypes. The truth variant sits on
    // haplotype 0, so hap 0 makes it a TP and hap 1 leaves the query call an FP -- both lanes
    // written, which is what distinguishes this from a single evaluate_variants call.
    GraphFixture f = build_fixture("ACGTACGT",
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}},
            {{2, 1, TYPE_SUB, "G", "A", GT_ALT1_REF, 60, 0, 0}});
    auto sc_data = make_superclusterData({"chr1"}, {8}, {2}, {f.sc}, f.ref);

    const std::vector< std::vector< std::vector<int> > > sc_groups = {{{0}, {0}}};
    precision_recall_wrapper(sc_data.get(), sc_groups, 0, 0, 1, false, false);

    EXPECT_EQ(ERRTYPE_TP, f.qvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_TP, f.tvars->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_FP, f.qvars->errtypes[HAP2][0]);
}

} // namespace
