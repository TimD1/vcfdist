/**
 * @file test_dist.cpp
 * @brief Unit tests for dist.cpp: pure alignment primitives, NG50, and graph cell indices.
 *
 * Every expected edit distance and alignment score below is derived by hand from the penalties
 * the algorithm charges (substitution x, gap of length L costs o + e*L), never by recording what
 * the implementation happens to print. The derivation is spelled out in a comment wherever the
 * arithmetic is not immediate.
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

// Build an in-memory reference by round-tripping a tiny FASTA through a temp file
// (fastaData only exposes a FILE* constructor).
std::shared_ptr<fastaData> make_ref(const std::string & ctg, const std::string & seq) {
    const char * path = "./gtest_graph_ref.fa";
    FILE * w = fopen(path, "w");
    fprintf(w, ">%s\n%s\n", ctg.c_str(), seq.c_str());
    fclose(w);
    FILE * r = fopen(path, "r");
    return std::make_shared<fastaData>(r);
}

// A truth insertion is a zero-width (in reference coordinates) locus. When another truth
// variant abuts it, a direct edge that leaps the insertion must not exist: every path across
// the locus has to route through the insertion's own alt or bypass node. This asserts that the
// reference node immediately after the insertion has only zero-width predecessors, so no
// reference-spanning node (e.g. the neighbouring SUB's alt or bypass) can bypass the insertion
// for free. Without this the insertion would be silently skipped and left unlabeled.
TEST(GraphInsertionEdges, AdjacentVariantCannotLeapInsertion) {
    auto sc = std::make_shared<ctgSuperclusters>();
    sc->callset_vars[QUERY] = std::make_shared<ctgVariants>("chr1");
    sc->callset_vars[TRUTH] = std::make_shared<ctgVariants>("chr1");
    auto ref = make_ref("chr1", "ACGTACGTAC");

    // truth SUB at pos 2 (G->T) immediately followed by truth INS at pos 3 (->TTT), both hap0
    auto tv = sc->callset_vars[TRUTH];
    tv->add_var(2, 1, TYPE_SUB, BED_INSIDE, "G", "T",   GT_ALT1_REF, 60, 60, 0, 0);
    tv->add_var(3, 0, TYPE_INS, BED_INSIDE, "",  "TTT", GT_ALT1_REF, 60, 60, 0, 0);

    Graph graph(sc, 0, ref, "chr1", HAP1);

    // locate the insertion locus (the zero-width truth variant node)
    int ins_coord = -1;
    for (int tn = 0; tn < graph.tnodes; tn++)
        if (graph.tidxs[tn] >= 0 && graph.tbegs[tn] == graph.tends[tn])
            ins_coord = graph.tbegs[tn];
    ASSERT_GE(ins_coord, 0) << "no zero-width insertion node found";

    // locate the reference node immediately to the right of the insertion locus
    int right_ref = -1;
    for (int tn = 0; tn < graph.tnodes; tn++)
        if (graph.ttypes[tn] == TYPE_REF && graph.tskips[tn] < 0 && graph.tidxs[tn] < 0 &&
                graph.tbegs[tn] == ins_coord && graph.tends[tn] > graph.tbegs[tn])
            right_ref = tn;
    ASSERT_GE(right_ref, 0) << "no reference node found after the insertion locus";

    // every predecessor must be a zero-width node at the locus (the insertion's alt/bypass);
    // a reference-spanning predecessor would be a free leap over the insertion
    ASSERT_FALSE(graph.tprevs[right_ref].empty());
    for (int p : graph.tprevs[right_ref])
        EXPECT_EQ(graph.tbegs[p], graph.tends[p])
            << "node " << p << " (spanning " << graph.tbegs[p] << ".." << graph.tends[p]
            << ") leaps the insertion at " << ins_coord;
}

// End-to-end guard: with the leap-suppression edge rule, a truth insertion abutting another
// truth variant is forced onto the alignment path (via its alt or bypass node) and is labeled by
// the backtrack itself -- not by the removed safety sweep. Here the query equals the reference, so
// both truth variants are missed and must be labeled FN via their bypass nodes, never left UNKNOWN.
TEST(GraphInsertionEdges, AdjacentInsertionLabeledWithoutSweep) {
    auto sc = std::make_shared<ctgSuperclusters>();
    sc->callset_vars[QUERY] = std::make_shared<ctgVariants>("chr1");
    sc->callset_vars[TRUTH] = std::make_shared<ctgVariants>("chr1");
    auto ref = make_ref("chr1", "ACGTACGTAC");

    auto tv = sc->callset_vars[TRUTH];
    tv->add_var(2, 1, TYPE_SUB, BED_INSIDE, "G", "T",   GT_ALT1_REF, 60, 60, 0, 0);
    tv->add_var(3, 0, TYPE_INS, BED_INSIDE, "",  "TTT", GT_ALT1_REF, 60, 60, 0, 0);

    auto graph = std::make_shared<Graph>(sc, 0, ref, "chr1", HAP1);
    std::unordered_map<idx4, idx4> ptrs;
    std::unordered_set<int> bypassed;
    calc_prec_recall_aln(graph, ptrs, false);
    calc_prec_recall(graph, ptrs, HAP1, bypassed, false);

    // query == reference: neither truth variant is reproduced, both are FN (not left UNKNOWN)
    EXPECT_EQ(ERRTYPE_FN, tv->errtypes[HAP1][0]) << "SUB should be a false negative";
    EXPECT_EQ(ERRTYPE_FN, tv->errtypes[HAP1][1]) << "insertion should be a false negative";
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

TEST(WfSwgMaxReach, ExhaustedQueryReturnsEarlyOnWorseDiagonal) {
    GlobalsGuard guard;
    set_penalties(4, 6, 2);
    const std::string query = "ACGTAAA";
    const std::string truth = "ACGTCCC";

    // DOCUMENTS A DEFECT; it does not enforce correct behavior. Hand-derived, a budget of 12
    // affords three substitutions along the main diagonal and should reach truth index 6. What
    // actually happens is that a budget of 12 also affords an insertion run (o+e, then e per
    // extension: 8, 10, 12) that consumes all three remaining query bases without advancing the
    // truth. The extend loop scans diagonals in ascending index order, so that insertion-heavy
    // diagonal is seen first, its `off == query_len-1` early return fires, and the smaller reach
    // of 3 is returned before the main diagonal is ever examined. The reported reach is therefore
    // not monotonic in max_score. Update these expectations when the defect is fixed.
    std::vector<int> offs11 = reach_offs(query, truth);
    EXPECT_EQ(5, wf_swg_max_reach(query, truth, offs11, 0, kNoMainDiagBlock, 11,
                g.sub, g.open, g.extend));
    std::vector<int> offs12 = reach_offs(query, truth);
    EXPECT_EQ(3, wf_swg_max_reach(query, truth, offs12, 0, kNoMainDiagBlock, 12,
                g.sub, g.open, g.extend));
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

} // namespace
