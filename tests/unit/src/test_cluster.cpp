/**
 * @file test_cluster.cpp
 * @brief Unit tests for cluster.cpp: index and range arithmetic, hap merging, and superclustering.
 */
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/cluster.h"
#include "../../../src/defs.h"
#include "../../../src/globals.h"
#include "test_helpers.h"

namespace {

/* Local helpers **********************************************************************************/

const int INT_MAXIMUM = std::numeric_limits<int>::max();

/**
 * @brief Builds a variant-free callset carrying the single trailing cluster boundary.
 *
 * The free functions in cluster.cpp read `clusters[cluster_start_indices[c]]` for every callset
 * even when that callset contributes no clusters, so an "empty" callset still needs one boundary.
 * @param[in] ctg Contig name
 * @return Variant container with n == 0, nc == 0, and clusters == {0}
 */
std::shared_ptr<ctgVariants> make_empty_callset(const std::string & ctg = "chr1") {
    std::shared_ptr<ctgVariants> vars = make_ctgVariants(ctg, {});
    set_clusters(vars, {0}, {}, {});
    return vars;
}

/* ctgSuperclusters ctor **************************************************************************/

TEST(CtgSuperclustersCtor, Size) {
    GlobalsGuard guard;
    ctgSuperclusters sc;
    EXPECT_EQ(size_t(CALLSETS), sc.callset_vars.size());
}

TEST(CtgSuperclustersCtor, Nulls) {
    GlobalsGuard guard;
    ctgSuperclusters sc;
    EXPECT_EQ(nullptr, sc.callset_vars[QUERY]);
    EXPECT_EQ(nullptr, sc.callset_vars[TRUTH]);
}

/* get_min_ref_pos ********************************************************************************/

TEST(GetMinRefPos, BothEmpty) {
    GlobalsGuard guard;
    std::shared_ptr<ctgSuperclusters> sc =
            make_ctgSuperclusters(make_empty_callset(), make_empty_callset());

    // both index ranges empty: min(INT_MAX, INT_MAX) - 1, asymmetric with get_max_ref_pos (#62)
    EXPECT_EQ(INT_MAXIMUM - 1, sc->get_min_ref_pos(0, 0, 0, 0));
}

TEST(GetMinRefPos, QueryOnly) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{40, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // one position left of the query variant
    EXPECT_EQ(39, sc->get_min_ref_pos(0, 1, 0, 0));
}

TEST(GetMinRefPos, TruthOnly) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{40, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(make_empty_callset(), tvars);
    EXPECT_EQ(39, sc->get_min_ref_pos(0, 0, 0, 1));
}

TEST(GetMinRefPos, BothQuerySmaller) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{20, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);
    EXPECT_EQ(9, sc->get_min_ref_pos(0, 1, 0, 1));
}

TEST(GetMinRefPos, BothTruthSmaller) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{20, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);
    EXPECT_EQ(9, sc->get_min_ref_pos(0, 1, 0, 1));
}

TEST(GetMinRefPos, UsesStartIndexOnly) {
    GlobalsGuard guard;

    // deliberately unsorted: a scan of [0, 2) would find position 5, an index read finds 100
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{100, 1, TYPE_SUB, "A", "C"}, {5, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // only poss[qvi_start] is read, so the sorted-input assumption is load-bearing
    EXPECT_EQ(99, sc->get_min_ref_pos(0, 2, 0, 0));
    EXPECT_NE(4, sc->get_min_ref_pos(0, 2, 0, 0));
}

TEST(GetMinRefPos, ContigStartClamped) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{0, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // the -1 flank is clamped away, since Graph() slices the reference from this offset (#167)
    EXPECT_EQ(0, sc->get_min_ref_pos(0, 1, 0, 0));
}

TEST(GetMinRefPos, ContigStartClampDoesNotHideSentinel) {
    GlobalsGuard guard;
    std::shared_ptr<ctgSuperclusters> sc =
            make_ctgSuperclusters(make_empty_callset(), make_empty_callset());

    // clamping the low end must leave the empty-range sentinel untouched
    EXPECT_EQ(INT_MAXIMUM - 1, sc->get_min_ref_pos(0, 0, 0, 0));
}

TEST(GetMinRefPos, Tie) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{30, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1",
            {{30, 5, TYPE_DEL, "ACGTA", "A"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);
    EXPECT_EQ(29, sc->get_min_ref_pos(0, 1, 0, 1));
}

/* get_max_ref_pos ********************************************************************************/

TEST(GetMaxRefPos, BothEmpty) {
    GlobalsGuard guard;
    std::shared_ptr<ctgSuperclusters> sc =
            make_ctgSuperclusters(make_empty_callset(), make_empty_callset());

    // the empty sentinel skips the trailing +1, unlike get_min_ref_pos's -1 (#62)
    EXPECT_EQ(INT_MAXIMUM, sc->get_max_ref_pos(0, 0, 0, 0));
    EXPECT_NE(INT_MAXIMUM - 1, sc->get_max_ref_pos(0, 0, 0, 0));
}

TEST(GetMaxRefPos, QuerySingle) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars =
            make_ctgVariants("chr1", {{40, 3, TYPE_DEL, "ACGT", "A"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // pos + rlen + 1
    EXPECT_EQ(44, sc->get_max_ref_pos(0, 1, 0, 0));
}

TEST(GetMaxRefPos, ScansFullRange) {
    GlobalsGuard guard;

    // the widest variant is in the middle of the range, so neither end index alone suffices
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C"}, {20, 30, TYPE_DEL, "A", "A"},
             {30, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // 20 + 30 + 1, not the first variant's 12 nor the last variant's 32
    EXPECT_EQ(51, sc->get_max_ref_pos(0, 3, 0, 0));
}

TEST(GetMaxRefPos, CrossCallset) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1",
            {{30, 5, TYPE_DEL, "ACGTA", "A"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // the truth variant reaches further right than the query variant
    EXPECT_EQ(36, sc->get_max_ref_pos(0, 1, 0, 1));
}

TEST(GetMaxRefPos, RlenZeroIns) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{7, 0, TYPE_INS, "", "ACGT"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // an insertion consumes no reference, so the end is one position right of its start
    EXPECT_EQ(8, sc->get_max_ref_pos(0, 1, 0, 0));
}

TEST(GetMaxRefPos, TruthOnly) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1",
            {{30, 5, TYPE_DEL, "ACGTA", "A"}, {50, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(make_empty_callset(), tvars);
    EXPECT_EQ(52, sc->get_max_ref_pos(0, 0, 0, 2));
}

/* get_supercluster_range *************************************************************************/

TEST(GetSuperclusterRange, BothEmptyDegenerate) {
    GlobalsGuard guard;
    std::shared_ptr<ctgSuperclusters> sc =
            make_ctgSuperclusters(make_empty_callset(), make_empty_callset());

    // neither callset contributes a cluster, so the initial values are returned unchanged
    std::vector<int> range = get_supercluster_range(sc->callset_vars, {0, 0}, {0, 0});
    ASSERT_EQ(size_t(2), range.size());
    EXPECT_EQ(INT_MAXIMUM, range[0]);
    EXPECT_EQ(-1, range[1]);
}

TEST(GetSuperclusterRange, QuerySingle) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1}, {5}, {15});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // one position left of pos, and one position right of pos+rlen
    std::vector<int> range = get_supercluster_range(sc->callset_vars, {0, 0}, {1, 0});
    EXPECT_EQ(9, range[0]);
    EXPECT_EQ(12, range[1]);
}

TEST(GetSuperclusterRange, MultiCluster) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C"}, {20, 1, TYPE_SUB, "A", "C"},
             {30, 1, TYPE_SUB, "A", "C"}, {40, 2, TYPE_DEL, "ACG", "A"}});
    set_clusters(qvars, {0, 2, 4}, {5, 25}, {25, 45});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // both clusters are included, so the range spans the first and last variants
    std::vector<int> range = get_supercluster_range(sc->callset_vars, {0, 0}, {2, 0});
    EXPECT_EQ(9, range[0]);
    EXPECT_EQ(43, range[1]);
}

TEST(GetSuperclusterRange, CrossCallsetMinMax) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{50, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1}, {45}, {55});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    set_clusters(tvars, {0, 1}, {5}, {15});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // truth supplies the minimum, query the maximum
    std::vector<int> range = get_supercluster_range(sc->callset_vars, {0, 0}, {1, 1});
    EXPECT_EQ(9, range[0]);
    EXPECT_EQ(52, range[1]);
}

TEST(GetSuperclusterRange, StartIdxOutOfBoundsErrors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C"}, {20, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1, 2}, {5, 15}, {15, 25});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // the guard is `start > nc`, so nc+1 == 3 is the first rejected start index
    EXPECT_EXIT(get_supercluster_range(sc->callset_vars, {3, 0}, {4, 0}),
            testing::ExitedWithCode(1), "Cluster start indices invalid");
}

TEST(GetSuperclusterRange, EndIdxOutOfBoundsErrors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C"}, {20, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1, 2}, {5, 15}, {15, 25});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // start == nc passes the `>` guard, so the end guard fires first; `>=` would invert this
    EXPECT_EXIT(get_supercluster_range(sc->callset_vars, {2, 0}, {3, 0}),
            testing::ExitedWithCode(1), "Cluster end indices invalid");
}

TEST(GetSuperclusterRange, OneCallsetEmpty) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1}, {5}, {15});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{1000, 1, TYPE_SUB, "A", "C"}});
    set_clusters(tvars, {0, 1}, {995}, {1005});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // truth has a variant but an empty cluster index range, so it is skipped entirely
    std::vector<int> range = get_supercluster_range(sc->callset_vars, {0, 0}, {1, 0});
    EXPECT_EQ(9, range[0]);
    EXPECT_EQ(12, range[1]);
}

/* get_next_variant_info **************************************************************************/

TEST(GetNextVariantInfo, BothExhausted) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{20, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // no callset has an unprocessed variant, so the result is flagged as not found
    var_info info = get_next_variant_info(sc->callset_vars, {{1, 1}}, {{1, 1}});
    EXPECT_FALSE(info.found);
    EXPECT_EQ(INT_MAXIMUM, info.start_pos);
    EXPECT_EQ(INT_MAXIMUM, info.end_pos);
}

TEST(GetNextVariantInfo, QueryOnly) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{5, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // the earlier truth variant is already consumed, so query wins despite being further right
    var_info info = get_next_variant_info(sc->callset_vars, {0, 1}, {1, 1});
    EXPECT_EQ(QUERY, info.callset_idx);
    EXPECT_EQ(10, info.start_pos);
    EXPECT_EQ(11, info.end_pos);
}

TEST(GetNextVariantInfo, TruthOnly) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{5, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1",
            {{30, 5, TYPE_DEL, "ACGTA", "A"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    var_info info = get_next_variant_info(sc->callset_vars, {1, 0}, {1, 1});
    EXPECT_EQ(TRUTH, info.callset_idx);
    EXPECT_EQ(30, info.start_pos);
    EXPECT_EQ(35, info.end_pos);
}

TEST(GetNextVariantInfo, QueryEarlier) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{20, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    var_info info = get_next_variant_info(sc->callset_vars, {0, 0}, {1, 1});
    EXPECT_EQ(QUERY, info.callset_idx);
    EXPECT_EQ(10, info.start_pos);
    EXPECT_EQ(11, info.end_pos);
}

TEST(GetNextVariantInfo, TruthEarlier) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{30, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{20, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    var_info info = get_next_variant_info(sc->callset_vars, {0, 0}, {1, 1});
    EXPECT_EQ(TRUTH, info.callset_idx);
    EXPECT_EQ(20, info.start_pos);
    EXPECT_EQ(21, info.end_pos);
}

TEST(GetNextVariantInfo, TiePrefersQuery) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{20, 3, TYPE_DEL, "ACGT", "A"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{20, 7, TYPE_DEL, "A", "A"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // QUERY is visited first and the comparison is strict `<`, so truth cannot displace it
    var_info info = get_next_variant_info(sc->callset_vars, {0, 0}, {1, 1});
    EXPECT_EQ(QUERY, info.callset_idx);
    EXPECT_EQ(20, info.start_pos);
    EXPECT_EQ(23, info.end_pos);
}

TEST(GetNextVariantInfo, EndPosOfWinner) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{12, 50, TYPE_DEL, "A", "A"}});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // the reported end position belongs to the leftmost-starting variant, not the widest one
    var_info info = get_next_variant_info(sc->callset_vars, {0, 0}, {1, 1});
    EXPECT_EQ(QUERY, info.callset_idx);
    EXPECT_EQ(11, info.end_pos);
    EXPECT_NE(62, info.end_pos);
}

/* get_supercluster_split_location ****************************************************************/

TEST(GetSuperclusterSplitLocation, FewerThanTwoEmpty) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", {{10, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1}, {5}, {15});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // a single variant cannot be split apart from anything
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {1, 0});
    EXPECT_TRUE(split.empty());
}

TEST(GetSuperclusterSplitLocation, NoGapsEmpty) {
    GlobalsGuard guard;

    // four abutting substitutions: every inter-variant gap is 0
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C"}, {11, 1, TYPE_SUB, "A", "C"},
             {12, 1, TYPE_SUB, "A", "C"}, {13, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1, 2, 3, 4}, {9, 10, 11, 12}, {11, 12, 13, 14});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // the score is gap/splits_to_halve_size, so a dense region scores 0 everywhere and never
    // beats the initial best of 0; an oversized dense supercluster is unsplittable (#63)
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {4, 0});
    EXPECT_TRUE(split.empty());
}

TEST(GetSuperclusterSplitLocation, SingleGap) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C"}, {11, 1, TYPE_SUB, "A", "C"},
             {50, 1, TYPE_SUB, "A", "C"}, {51, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 1, 2, 3, 4}, {9, 10, 49, 50}, {11, 12, 51, 52});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // the only positive gap is 12..50, which happens to split the range exactly in half
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {4, 0});
    ASSERT_EQ(size_t(CALLSETS), split.size());
    EXPECT_EQ(2, split[idx(QUERY)]);
    EXPECT_EQ(0, split[idx(TRUTH)]);
}

TEST(GetSuperclusterSplitLocation, PrefersCentralGap) {
    GlobalsGuard guard;

    // two 20bp gaps, 0..20 and 20..40; the supercluster range is [-1, 42)
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{0, 0, TYPE_INS, "", "A"}, {20, 0, TYPE_INS, "", "A"},
             {40, 0, TYPE_INS, "", "A"}, {41, 0, TYPE_INS, "", "A"}});
    set_clusters(qvars, {0, 1, 2, 3, 4}, {-1, 19, 39, 40}, {1, 21, 41, 42});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // equal gaps, so the tie is broken by size_reduction_factor: the second gap's midpoint 30
    // is nearer the range centre 20.5 than the first gap's midpoint 10
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {4, 0});
    ASSERT_EQ(size_t(CALLSETS), split.size());
    EXPECT_EQ(2, split[idx(QUERY)]);
}

TEST(GetSuperclusterSplitLocation, EdgeGapLowScore) {
    GlobalsGuard guard;

    // a 30bp gap at the left edge, then two 10bp gaps; the last variant is 100bp wide, so the
    // supercluster range is [-1, 151) and the left gap is heavily off-centre
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{0, 0, TYPE_INS, "", "A"}, {30, 0, TYPE_INS, "", "A"},
             {40, 0, TYPE_INS, "", "A"}, {50, 100, TYPE_DEL, "A", "A"}});
    set_clusters(qvars, {0, 1, 2, 3, 4}, {-1, 29, 39, 49}, {1, 31, 41, 151});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // the 30bp edge gap scores ~4.81 while the central 10bp gap before variant 3 scores ~5.20,
    // so the larger but off-centre gap loses
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {4, 0});
    ASSERT_EQ(size_t(CALLSETS), split.size());
    EXPECT_EQ(3, split[idx(QUERY)]);
    EXPECT_NE(1, split[idx(QUERY)]);
}

TEST(GetSuperclusterSplitLocation, Log2DivZeroGuard) {
    GlobalsGuard guard;

    // the first variant is 30bp wide and ends at 30, exactly the range end; the second variant is
    // nested inside it, so the gap is 0 and the candidate split midpoint lands on the range end
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{0, 30, TYPE_DEL, "A", "A"}, {25, 4, TYPE_DEL, "ACGTA", "A"}});
    set_clusters(qvars, {0, 1, 2}, {-1, 24}, {30, 30});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // size_reduction_factor is exactly 31/31 == 1.0, whose log2 is 0; the `< 1.0` and
    // `log2_srf != 0` guards keep -1/log2_srf from dividing by zero, so the score stays 0 and
    // no split location is reported (#63)
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {2, 0});
    EXPECT_TRUE(split.empty());
}

TEST(GetSuperclusterSplitLocation, IndexAfterCurr) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{0, 0, TYPE_INS, "", "A"}, {10, 0, TYPE_INS, "", "A"},
             {100, 0, TYPE_INS, "", "A"}});
    set_clusters(qvars, {0, 1, 2, 3}, {-1, 9, 99}, {1, 11, 101});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());

    // the winning gap is 10..100, between variants 1 and 2; the returned index is that of the
    // first variant to the RIGHT of the split, because split_indices is incremented past curr
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {3, 0});
    ASSERT_EQ(size_t(CALLSETS), split.size());
    EXPECT_EQ(2, split[idx(QUERY)]);
    EXPECT_NE(1, split[idx(QUERY)]);
}

TEST(GetSuperclusterSplitLocation, CrossCallsetGap) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{0, 0, TYPE_INS, "", "A"}, {50, 0, TYPE_INS, "", "A"}});
    set_clusters(qvars, {0, 1, 2}, {-1, 49}, {1, 51});
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {{10, 0, TYPE_INS, "", "A"}});
    set_clusters(tvars, {0, 1}, {9}, {11});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);

    // the winning 40bp gap runs from the truth variant at 10 to the query variant at 50, so the
    // split index has advanced on both callsets
    std::vector<int> split = get_supercluster_split_location(sc->callset_vars, {0, 0}, {2, 1});
    ASSERT_EQ(size_t(CALLSETS), split.size());
    EXPECT_EQ(1, split[idx(QUERY)]);
    EXPECT_EQ(1, split[idx(TRUTH)]);
}

/* split_cluster **********************************************************************************/

/**
 * @brief Builds the six-variant query callset shared by the split_cluster tests.
 *
 * Variants sit at positions 10 through 60 in steps of 10, grouped into three clusters of two.
 * @return Variant container with n == 6, nc == 3, and clusters == {0, 2, 4, 6}
 */
std::shared_ptr<ctgVariants> make_split_cluster_query() {
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C"}, {20, 1, TYPE_SUB, "A", "C"},
             {30, 1, TYPE_SUB, "A", "C"}, {40, 1, TYPE_SUB, "A", "C"},
             {50, 1, TYPE_SUB, "A", "C"}, {60, 1, TYPE_SUB, "A", "C"}});
    set_clusters(qvars, {0, 2, 4, 6}, {5, 25, 45}, {25, 45, 65});
    return qvars;
}

TEST(SplitCluster, AtExistingBoundaryNoop) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1",
            {{100, 1, TYPE_SUB, "A", "C"}, {110, 1, TYPE_SUB, "A", "C"}});
    set_clusters(tvars, {0, 1, 2}, {95, 105}, {105, 115});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints = {{0, 0}, {3, 2}};

    // query variant 2 and truth variant 1 are both already cluster boundaries
    EnumArray<callset_t, int, CALLSET_SLOTS> split = split_cluster(sc->callset_vars, {2, 1}, breakpoints, 0);

    EXPECT_EQ(1, split[QUERY]);
    EXPECT_EQ(1, split[TRUTH]);
    EXPECT_EQ(3, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 4, 6}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({5, 25, 45}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({25, 45, 65}), qvars->right_reaches);
    EXPECT_EQ(2, tvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2}), tvars->clusters);
    EXPECT_EQ(std::vector<int>({105, 115}), tvars->right_reaches);

    // nothing was inserted, so no later breakpoint moved
    EXPECT_EQ(0, breakpoints[0][QUERY]);
    EXPECT_EQ(0, breakpoints[0][TRUTH]);
    EXPECT_EQ(3, breakpoints[1][QUERY]);
    EXPECT_EQ(2, breakpoints[1][TRUTH]);
}

TEST(SplitCluster, MidclusterInserts) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints = {{0, 0}, {3, 0}};

    // query variant 3 sits inside cluster 1 (variants 2..3), so that cluster is split in two
    EnumArray<callset_t, int, CALLSET_SLOTS> split = split_cluster(sc->callset_vars, {3, 0}, breakpoints, 0);

    EXPECT_EQ(2, split[QUERY]);
    EXPECT_EQ(0, split[TRUTH]);
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 3, 4, 6}), qvars->clusters);
    EXPECT_EQ(size_t(4), qvars->left_reaches.size());
    EXPECT_EQ(size_t(4), qvars->right_reaches.size());

    // the truth callset holds no variants, so its lanes are untouched
    EXPECT_EQ(0, sc->callset_vars[TRUTH]->nc);
    EXPECT_EQ(std::vector<int>({0}), sc->callset_vars[TRUTH]->clusters);
}

TEST(SplitCluster, ReachReassignment) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints = {{0, 0}, {3, 0}};

    split_cluster(sc->callset_vars, {3, 0}, breakpoints, 0);

    // the left piece keeps its left reach but its right reach becomes the split position, while
    // the new right piece starts at the split position and inherits the original right reach
    EXPECT_EQ(std::vector<int>({5, 25, 40, 45}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({25, 40, 45, 65}), qvars->right_reaches);
}

TEST(SplitCluster, LaterBreakpointsIncremented) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints = {{0, 0}, {2, 0}, {3, 0}};

    split_cluster(sc->callset_vars, {3, 0}, breakpoints, 0);

    // inserting a cluster shifts every later cluster index on the split callset by one
    EXPECT_EQ(0, breakpoints[0][QUERY]);
    EXPECT_EQ(0, breakpoints[0][TRUTH]);
    EXPECT_EQ(3, breakpoints[1][QUERY]);
    EXPECT_EQ(0, breakpoints[1][TRUTH]);
    EXPECT_EQ(4, breakpoints[2][QUERY]);
    EXPECT_EQ(0, breakpoints[2][TRUTH]);

    // only breakpoints strictly after breakpoint_idx move
    std::shared_ptr<ctgVariants> qvars2 = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc2 = make_ctgSuperclusters(qvars2, make_empty_callset());
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints2 = {{0, 0}, {2, 0}, {3, 0}};
    split_cluster(sc2->callset_vars, {3, 0}, breakpoints2, 1);
    EXPECT_EQ(2, breakpoints2[1][QUERY]);
    EXPECT_EQ(0, breakpoints2[1][TRUTH]);
    EXPECT_EQ(4, breakpoints2[2][QUERY]);
    EXPECT_EQ(0, breakpoints2[2][TRUTH]);
}

TEST(SplitCluster, BothCallsetsIndependent) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1",
            {{100, 1, TYPE_SUB, "A", "C"}, {110, 1, TYPE_SUB, "A", "C"}});
    set_clusters(tvars, {0, 2}, {95}, {115});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints = {{0, 0}, {3, 1}};

    // both callsets split mid-cluster, at different cluster indices
    EnumArray<callset_t, int, CALLSET_SLOTS> split = split_cluster(sc->callset_vars, {3, 1}, breakpoints, 0);

    EXPECT_EQ(2, split[QUERY]);
    EXPECT_EQ(1, split[TRUTH]);
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 3, 4, 6}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({5, 25, 40, 45}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({25, 40, 45, 65}), qvars->right_reaches);
    EXPECT_EQ(2, tvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2}), tvars->clusters);
    EXPECT_EQ(std::vector<int>({95, 110}), tvars->left_reaches);
    EXPECT_EQ(std::vector<int>({110, 115}), tvars->right_reaches);

    // each callset's later breakpoints are incremented independently
    EXPECT_EQ(4, breakpoints[1][QUERY]);
    EXPECT_EQ(2, breakpoints[1][TRUTH]);
}

TEST(SplitCluster, ClustIdxMinusOneSafety) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints = {{0, 0}, {3, 0}};

    // variant 0 equals clusters[0], so the equality branch is taken and right_reaches[clust_idx-1]
    // is never indexed at -1; callers only pass indices at or after the first boundary
    EnumArray<callset_t, int, CALLSET_SLOTS> split = split_cluster(sc->callset_vars, {0, 0}, breakpoints, 0);

    EXPECT_EQ(0, split[QUERY]);
    EXPECT_EQ(0, split[TRUTH]);
    EXPECT_EQ(3, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 4, 6}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({5, 25, 45}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({25, 45, 65}), qvars->right_reaches);
    EXPECT_EQ(3, breakpoints[1][QUERY]);
    EXPECT_EQ(0, breakpoints[1][TRUTH]);
}

/* split_large_supercluster ***********************************************************************/

/**
 * @brief Builds a four-variant query callset with one cluster per variant.
 *
 * Variants sit at positions 0, 100, 200 and 300, so the supercluster spans 302+rlen bases.
 * @param[in] rlen Reference length of every variant
 * @return Variant container with n == 4, nc == 4, and clusters == {0, 1, 2, 3, 4}
 */
std::shared_ptr<ctgVariants> make_spread_query(int rlen) {
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{0, rlen, TYPE_DEL, "A", "A"}, {100, rlen, TYPE_DEL, "A", "A"},
             {200, rlen, TYPE_DEL, "A", "A"}, {300, rlen, TYPE_DEL, "A", "A"}});
    set_clusters(qvars, {0, 1, 2, 3, 4}, {-1, 99, 199, 299},
            {rlen + 1, rlen + 101, rlen + 201, rlen + 301});
    return qvars;
}

/**
 * @brief Builds a four-variant query callset held in a single cluster.
 *
 * Splitting this supercluster must insert a new cluster boundary, unlike make_spread_query().
 * @return Variant container with n == 4, nc == 1, and clusters == {0, 4}
 */
std::shared_ptr<ctgVariants> make_one_cluster_query() {
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1",
            {{0, 0, TYPE_INS, "", "A"}, {100, 0, TYPE_INS, "", "A"},
             {200, 0, TYPE_INS, "", "A"}, {300, 0, TYPE_INS, "", "A"}});
    set_clusters(qvars, {0, 4}, {-1}, {301});
    return qvars;
}

TEST(SplitLargeSupercluster, AlreadySmallNoop) {
    GlobalsGuard guard;
    g.max_supercluster_size = 1000;
    std::shared_ptr<ctgVariants> qvars = make_spread_query(0);
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    EnumArray<callset_t, int, CALLSET_SLOTS> end_indices = {{4, 0}};

    // the supercluster spans 302bp, well within the limit
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    ASSERT_EQ(size_t(2), breakpoints.size());
    EXPECT_EQ(0, breakpoints[0][QUERY]);
    EXPECT_EQ(0, breakpoints[0][TRUTH]);
    EXPECT_EQ(4, breakpoints[1][QUERY]);
    EXPECT_EQ(0, breakpoints[1][TRUTH]);
    EXPECT_EQ(4, end_indices[QUERY]);
    EXPECT_EQ(0, end_indices[TRUTH]);
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2, 3, 4}), qvars->clusters);
}

TEST(SplitLargeSupercluster, OneSplit) {
    GlobalsGuard guard;
    g.max_supercluster_size = 150;
    std::shared_ptr<ctgVariants> qvars = make_spread_query(0);
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    EnumArray<callset_t, int, CALLSET_SLOTS> end_indices = {{4, 0}};

    // 302bp exceeds the limit, and the central gap splits it into two 102bp pieces
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    ASSERT_EQ(size_t(3), breakpoints.size());
    EXPECT_EQ(0, breakpoints[0][QUERY]);
    EXPECT_EQ(0, breakpoints[0][TRUTH]);
    EXPECT_EQ(2, breakpoints[1][QUERY]);
    EXPECT_EQ(0, breakpoints[1][TRUTH]);
    EXPECT_EQ(4, breakpoints[2][QUERY]);
    EXPECT_EQ(0, breakpoints[2][TRUTH]);

    // the split fell on an existing cluster boundary, so no cluster was inserted
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2, 3, 4}), qvars->clusters);
    EXPECT_EQ(4, end_indices[QUERY]);
    EXPECT_EQ(0, end_indices[TRUTH]);
}

TEST(SplitLargeSupercluster, MultipleSplits) {
    GlobalsGuard guard;
    g.max_supercluster_size = 80;
    std::shared_ptr<ctgVariants> qvars = make_spread_query(0);
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    EnumArray<callset_t, int, CALLSET_SLOTS> end_indices = {{4, 0}};

    // each 102bp half still exceeds 80bp, so a second round splits both of them
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    ASSERT_EQ(size_t(5), breakpoints.size());
    for (size_t i = 0; i < breakpoints.size(); i++) {
        EXPECT_EQ(int(i), breakpoints[i][QUERY]) << "breakpoint " << i;
        EXPECT_EQ(0, breakpoints[i][TRUTH]) << "breakpoint " << i;
    }

    // breakpoints are strictly increasing on the query callset
    for (size_t i = 1; i < breakpoints.size(); i++) {
        EXPECT_LT(breakpoints[i-1][QUERY], breakpoints[i][QUERY]) << "breakpoint " << i;
    }
}

TEST(SplitLargeSupercluster, NoValidSplitBails) {
    GlobalsGuard guard;
    g.max_supercluster_size = 150;

    // four 100bp deletions abutting end to end: 402bp wide, with no gap anywhere to split at
    std::shared_ptr<ctgVariants> qvars = make_spread_query(100);
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    EnumArray<callset_t, int, CALLSET_SLOTS> end_indices = {{4, 0}};

    testing::internal::CaptureStderr();
    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);
    std::string err = testing::internal::GetCapturedStderr();

    // the oversized supercluster is retained unchanged
    ASSERT_EQ(size_t(2), breakpoints.size());
    EXPECT_EQ(0, breakpoints[0][QUERY]);
    EXPECT_EQ(0, breakpoints[0][TRUTH]);
    EXPECT_EQ(4, breakpoints[1][QUERY]);
    EXPECT_EQ(0, breakpoints[1][TRUTH]);
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2, 3, 4}), qvars->clusters);

    // and the condition is surfaced rather than passed over silently (#64)
    EXPECT_NE(std::string::npos, err.find("No valid split location for oversized supercluster"));
    EXPECT_NE(std::string::npos, err.find("size 402"));
    EXPECT_NE(std::string::npos, err.find("--max-supercluster-size 150"));
    EXPECT_NE(std::string::npos, err.find("retaining"));
}

TEST(SplitLargeSupercluster, MutatesEndIndices) {
    GlobalsGuard guard;
    g.max_supercluster_size = 150;
    std::shared_ptr<ctgVariants> qvars = make_one_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    EnumArray<callset_t, int, CALLSET_SLOTS> end_indices = {{1, 0}};

    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    // the in/out end indices are overwritten with the final breakpoint, which the insertion moved
    EXPECT_EQ(2, end_indices[QUERY]);
    EXPECT_EQ(0, end_indices[TRUTH]);
    const auto & last = breakpoints[breakpoints.size()-1];
    EXPECT_EQ(last[QUERY], end_indices[QUERY]);
    EXPECT_EQ(last[TRUTH], end_indices[TRUTH]);
}

TEST(SplitLargeSupercluster, BreakpointShiftAfterInsert) {
    GlobalsGuard guard;
    g.max_supercluster_size = 150;
    std::shared_ptr<ctgVariants> qvars = make_one_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    EnumArray<callset_t, int, CALLSET_SLOTS> end_indices = {{1, 0}};

    std::vector< EnumArray<callset_t, int, CALLSET_SLOTS> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    // splitting at variant 2 inserted a cluster, so the trailing breakpoint shifted 1 -> 2
    ASSERT_EQ(size_t(3), breakpoints.size());
    EXPECT_EQ(0, breakpoints[0][QUERY]);
    EXPECT_EQ(0, breakpoints[0][TRUTH]);
    EXPECT_EQ(1, breakpoints[1][QUERY]);
    EXPECT_EQ(0, breakpoints[1][TRUTH]);
    EXPECT_EQ(2, breakpoints[2][QUERY]);
    EXPECT_EQ(0, breakpoints[2][TRUTH]);

    // the full post-state of the mutated cluster lanes
    EXPECT_EQ(2, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 4}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({-1, 200}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({200, 301}), qvars->right_reaches);
}

/* load_and_merge_callset_vars_across_haps ********************************************************/

/**
 * @brief Builds a one-contig superclusterData ready to receive a merged callset.
 *
 * The merge writes through `superclusters[ctg]`, so that entry has to exist before it runs.
 * @param[in] ctg Contig name
 * @param[in] length Contig length
 * @return Supercluster data holding one contig with both callsets empty
 */
std::shared_ptr<superclusterData> make_merge_target(const std::string & ctg = "chr1",
        int length = 1000) {
    std::shared_ptr<ctgSuperclusters> sc =
            make_ctgSuperclusters(make_empty_callset(ctg), make_empty_callset(ctg));
    return make_superclusterData({ctg}, {length}, {sc});
}

/**
 * @brief Wraps two per-haplotype variant containers in the contig-map form the merge expects.
 * @param[in] hap1 HAP1 variants
 * @param[in] hap2 HAP2 variants
 * @param[in] ctg Contig name
 * @return Per-haplotype maps from contig name to variants
 */
EnumArray<hap_t, std::unordered_map< std::string, std::shared_ptr<ctgVariants> >, HAP_SLOTS> make_hap_vars(
        std::shared_ptr<ctgVariants> hap1, std::shared_ptr<ctgVariants> hap2,
        const std::string & ctg = "chr1") {
    EnumArray<hap_t, std::unordered_map< std::string, std::shared_ptr<ctgVariants> >, HAP_SLOTS> vars;
    vars[HAP1][ctg] = hap1;
    vars[HAP2][ctg] = hap2;
    return vars;
}

/**
 * @brief Returns the merged QUERY callset written onto a contig by the merge.
 * @param[in] sc_data Supercluster data the merge wrote into
 * @param[in] ctg Contig name
 * @return Merged query variants for that contig
 */
std::shared_ptr<ctgVariants> merged_query(std::shared_ptr<superclusterData> sc_data,
        const std::string & ctg = "chr1") {
    return sc_data->superclusters[ctg]->callset_vars[QUERY];
}

/**
 * @brief Asserts that every cluster ahead of the trailing sentinel holds at least one variant.
 * @param[in] merged Merged callset written by the merge, on a non-empty contig
 */
void expect_no_empty_clusters(std::shared_ptr<ctgVariants> merged) {
    ASSERT_EQ(size_t(merged->nc), merged->clusters.size());
    ASSERT_GE(merged->nc, 2);
    EXPECT_EQ(merged->n, merged->clusters.back()) << "sentinel must start at n";
    for (int c = 0; c + 1 < merged->nc; c++) {
        EXPECT_LT(merged->clusters[c], merged->clusters[c+1]) << "cluster " << c << " is empty";
    }
}

TEST(LoadAndMerge, EmptyContigKeepsTrailingBoundary) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();
    auto vars = make_hap_vars(make_ctgVariants("chr1", {}), make_ctgVariants("chr1", {}));

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the early `continue` stores a container holding no cluster, but still carrying the single
    // trailing boundary that supercluster() reads for every callset; leaving clusters empty here
    // segfaulted on any contig the other callset had variants on (#166)
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    EXPECT_EQ(0, merged->n);
    EXPECT_EQ(0, merged->nc);
    EXPECT_EQ(std::vector<int>({0}), merged->clusters);
    EXPECT_EQ(merged->n, merged->clusters.back()) << "boundary must start at n";

    // nc stays 0, as wf_swg_cluster() also leaves it, so clusters holds nc+1 entries here while a
    // non-empty contig counts its sentinel and holds nc; the reaches are never indexed either way
    EXPECT_EQ(size_t(merged->nc) + 1, merged->clusters.size());
    EXPECT_TRUE(merged->left_reaches.empty());
    EXPECT_TRUE(merged->right_reaches.empty());
}

TEST(LoadAndMerge, HomVariant) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();

    // identical position, REF and ALT on both haplotypes, but opposite single-hap genotypes
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {5}, {15});
    std::shared_ptr<ctgVariants> hap2 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_REF_ALT}});
    set_clusters(hap2, {0, 1}, {5}, {15});
    auto vars = make_hap_vars(hap1, hap2);

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the two are collapsed into one homozygous record, and both haps advance together
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(1, merged->n);
    EXPECT_EQ(10, merged->poss[0]);
    EXPECT_EQ(GT_ALT_ALT, merged->orig_gts[0]);
}

TEST(LoadAndMerge, HetHap1First) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {5}, {12});
    std::shared_ptr<ctgVariants> hap2 = make_ctgVariants("chr1",
            {{20, 1, TYPE_SUB, "A", "G", GT_REF_ALT}});
    set_clusters(hap2, {0, 1}, {15}, {25});
    auto vars = make_hap_vars(hap1, hap2);

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the leftmost variant is emitted first and each keeps its own single-hap genotype
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(2, merged->n);
    EXPECT_EQ(std::vector<int>({10, 20}), merged->poss);
    EXPECT_EQ(GT_ALT_REF, merged->orig_gts[0]);
    EXPECT_EQ(GT_REF_ALT, merged->orig_gts[1]);

    // the reaches do not touch, so the two clusters stay separate ahead of the sentinel
    EXPECT_EQ(std::vector<int>({0, 1, 2}), merged->clusters);
    EXPECT_EQ(std::vector<int>({5, 15, INT_MAXIMUM}), merged->left_reaches);
    EXPECT_EQ(std::vector<int>({12, 25, INT_MAXIMUM}), merged->right_reaches);
}

TEST(LoadAndMerge, HetTiePrefersIns) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();

    // both haplotypes carry a variant at position 20, and HAP2's is the insertion
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{20, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {15}, {25});
    std::shared_ptr<ctgVariants> hap2 = make_ctgVariants("chr1",
            {{20, 0, TYPE_INS, "", "GG", GT_REF_ALT}});
    set_clusters(hap2, {0, 1}, {15}, {25});
    auto vars = make_hap_vars(hap1, hap2);

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // an insertion consumes no reference, so it is emitted ahead of the co-located substitution
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(2, merged->n);
    EXPECT_EQ(TYPE_INS, merged->types[0]);
    EXPECT_EQ("GG", merged->alts[0]);
    EXPECT_EQ(TYPE_SUB, merged->types[1]);
}

TEST(LoadAndMerge, HetTieDefaultHap1) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();

    // co-located but differing ALTs, and neither is an insertion
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{20, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {15}, {25});
    std::shared_ptr<ctgVariants> hap2 = make_ctgVariants("chr1",
            {{20, 1, TYPE_SUB, "A", "G", GT_REF_ALT}});
    set_clusters(hap2, {0, 1}, {15}, {25});
    auto vars = make_hap_vars(hap1, hap2);

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the INS tie-break does not apply and HAP2 is not strictly left of HAP1, so HAP1 wins
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(2, merged->n);
    EXPECT_EQ("C", merged->alts[0]);
    EXPECT_EQ("G", merged->alts[1]);
}

TEST(LoadAndMerge, Hap1Only) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF},
             {12, 1, TYPE_SUB, "A", "G", GT_ALT_REF}});
    set_clusters(hap1, {0, 2}, {5}, {17});
    auto vars = make_hap_vars(hap1, make_ctgVariants("chr1", {}));

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // HAP2 contributes nothing, so the HAP1 cluster passes through unchanged
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(2, merged->n);
    EXPECT_EQ(std::vector<int>({10, 12}), merged->poss);
    EXPECT_EQ(std::vector<int>({0, 2}), merged->clusters);
    EXPECT_EQ(std::vector<int>({5, INT_MAXIMUM}), merged->left_reaches);
    EXPECT_EQ(std::vector<int>({17, INT_MAXIMUM}), merged->right_reaches);
}

TEST(LoadAndMerge, Hap2Only) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();
    std::shared_ptr<ctgVariants> hap2 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_REF_ALT},
             {12, 1, TYPE_SUB, "A", "G", GT_REF_ALT}});
    set_clusters(hap2, {0, 2}, {5}, {17});
    auto vars = make_hap_vars(make_ctgVariants("chr1", {}), hap2);

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the HAP2-only initialization branch mirrors the HAP1-only one
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(2, merged->n);
    EXPECT_EQ(std::vector<int>({10, 12}), merged->poss);
    EXPECT_EQ(std::vector<int>({0, 2}), merged->clusters);
    EXPECT_EQ(std::vector<int>({5, INT_MAXIMUM}), merged->left_reaches);
    EXPECT_EQ(std::vector<int>({17, INT_MAXIMUM}), merged->right_reaches);
}

TEST(LoadAndMerge, TwoClustersSeparate) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF},
             {110, 1, TYPE_SUB, "A", "G", GT_ALT_REF}});
    set_clusters(hap1, {0, 1, 2}, {5, 105}, {15, 115});
    auto vars = make_hap_vars(hap1, make_ctgVariants("chr1", {}));

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // 15 < 105, so the second cluster starts a new supercluster rather than extending the first
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    EXPECT_EQ(std::vector<int>({0, 1, 2}), merged->clusters);
    EXPECT_EQ(std::vector<int>({5, 105, INT_MAXIMUM}), merged->left_reaches);
    EXPECT_EQ(std::vector<int>({15, 115, INT_MAXIMUM}), merged->right_reaches);
}

TEST(LoadAndMerge, TwoClustersOverlap) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();

    // HAP2's cluster reaches left to 12, inside HAP1's cluster reaching right to 20
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {5}, {20});
    std::shared_ptr<ctgVariants> hap2 = make_ctgVariants("chr1",
            {{14, 1, TYPE_SUB, "A", "G", GT_REF_ALT}});
    set_clusters(hap2, {0, 1}, {12}, {25});
    auto vars = make_hap_vars(hap1, hap2);

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the overlapping clusters collapse into one spanning the union of both reaches
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(2, merged->n);
    EXPECT_EQ(std::vector<int>({0, 2}), merged->clusters);
    EXPECT_EQ(std::vector<int>({5, INT_MAXIMUM}), merged->left_reaches);
    EXPECT_EQ(std::vector<int>({25, INT_MAXIMUM}), merged->right_reaches);
}

TEST(LoadAndMerge, SentinelAppended) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {5}, {15});
    auto vars = make_hap_vars(hap1, make_ctgVariants("chr1", {}));

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the trailing sentinel starts at n and reaches int::max in both directions, and nc counts it
    // (unlike wf_swg_cluster's nc, which does not) -- so clusters.size() == nc here
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    EXPECT_EQ(size_t(merged->nc), merged->clusters.size());
    EXPECT_EQ(merged->n, merged->clusters.back());
    EXPECT_EQ(INT_MAXIMUM, merged->left_reaches.back());
    EXPECT_EQ(INT_MAXIMUM, merged->right_reaches.back());

    // the last variant closes its cluster with no cluster left to start, so curr_* is left holding
    // an empty (int::max, int::min) pair that must not be saved ahead of the sentinel (#168)
    EXPECT_EQ(2, merged->nc);
    EXPECT_EQ(std::vector<int>({0, 1}), merged->clusters);
    EXPECT_EQ(std::vector<int>({5, INT_MAXIMUM}), merged->left_reaches);
    EXPECT_EQ(std::vector<int>({15, INT_MAXIMUM}), merged->right_reaches);
    expect_no_empty_clusters(merged);
}

TEST(LoadAndMerge, SingleVariantContigHomBothHaps) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();

    // as HomVariant, but asserting the cluster lanes: both haplotypes exhaust on the same
    // iteration here, rather than one running out ahead of the other
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {5}, {15});
    std::shared_ptr<ctgVariants> hap2 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_REF_ALT}});
    set_clusters(hap2, {0, 1}, {5}, {15});
    auto vars = make_hap_vars(hap1, hap2);

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);

    // the collapsed homozygous record is the contig's only variant, so this is the same
    // single-cluster shape as SentinelAppended (#168)
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(1, merged->n);
    EXPECT_EQ(2, merged->nc);
    EXPECT_EQ(std::vector<int>({0, 1}), merged->clusters);
    EXPECT_EQ(std::vector<int>({5, INT_MAXIMUM}), merged->left_reaches);
    EXPECT_EQ(std::vector<int>({15, INT_MAXIMUM}), merged->right_reaches);
    expect_no_empty_clusters(merged);
}

TEST(LoadAndMerge, SingleVariantContigSuperclustersCorrectly) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_merge_target();
    std::shared_ptr<ctgVariants> hap1 = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF}});
    set_clusters(hap1, {0, 1}, {5}, {15});
    auto vars = make_hap_vars(hap1, make_ctgVariants("chr1", {}));

    sc_data->load_and_merge_callset_vars_across_haps(QUERY, vars);
    sc_data->supercluster(/* print = */ false);

    // the sentinel's int::max left reach is what ends the supercluster scan, so dropping the empty
    // cluster leaves the one variant in supercluster 0 exactly as before
    std::shared_ptr<ctgVariants> merged = merged_query(sc_data);
    ASSERT_EQ(size_t(1), merged->superclusters.size());
    EXPECT_EQ(0, merged->superclusters[0]);
}

/* supercluster ***********************************************************************************/

/**
 * @brief Builds a callset in the merged form supercluster() consumes.
 *
 * load_and_merge_callset_vars_across_haps() counts its trailing sentinel in nc, so all three
 * cluster lanes hold nc entries -- unlike wf_swg_cluster's output, where clusters holds nc+1 and
 * the reaches hold nc. The sentinel is appended here so callers describe only real clusters.
 * @param[in] vars Variants to hold, in ascending position order
 * @param[in] clusters Variant index at the start of each real cluster
 * @param[in] left_reaches Leftmost reach of each real cluster
 * @param[in] right_reaches Rightmost reach of each real cluster
 * @param[in] ctg Contig name
 * @return Variant container carrying a sentinel-terminated clustering
 */
std::shared_ptr<ctgVariants> make_merged_callset(const std::vector<var_desc> & vars,
        const std::vector<int> & clusters, const std::vector<int> & left_reaches,
        const std::vector<int> & right_reaches, const std::string & ctg = "chr1") {
    std::shared_ptr<ctgVariants> merged = make_ctgVariants(ctg, vars);
    std::vector<int> cl = clusters, lr = left_reaches, rr = right_reaches;
    cl.push_back(merged->n);
    lr.push_back(INT_MAXIMUM);
    rr.push_back(INT_MAXIMUM);
    set_clusters(merged, cl, lr, rr, int(cl.size()));
    return merged;
}

/**
 * @brief Wraps one contig's callsets in a superclusterData and superclusters them.
 * @param[in,out] qvars Query variants, whose supercluster lane is filled in
 * @param[in,out] tvars Truth variants, whose supercluster lane is filled in
 * @param[in] ctg Contig name
 * @return The superclustered data, kept alive for the caller's assertions
 * @throws WARNING if a supercluster exceeds g.max_supercluster_size
 */
std::shared_ptr<superclusterData> run_supercluster(std::shared_ptr<ctgVariants> qvars,
        std::shared_ptr<ctgVariants> tvars, const std::string & ctg = "chr1") {
    std::shared_ptr<superclusterData> sc_data =
            make_superclusterData({ctg}, {1000}, {make_ctgSuperclusters(qvars, tvars)});
    sc_data->supercluster();
    return sc_data;
}

/**
 * @brief Describes one single-base substitution at each of the given positions.
 * @param[in] poss Reference positions, ascending
 * @return Variant descriptors ready for make_ctgVariants()
 */
std::vector<var_desc> subs_at(const std::vector<int> & poss) {
    std::vector<var_desc> vars;
    for (int pos : poss) vars.push_back({pos, 1, TYPE_SUB, "A", "C"});
    return vars;
}

TEST(Supercluster, EmptySkipped) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_empty_callset();
    std::shared_ptr<ctgVariants> tvars = make_empty_callset();

    // no variants on either callset, so the contig is skipped before any cluster is read
    run_supercluster(qvars, tvars);

    EXPECT_EQ(0, qvars->n);
    EXPECT_TRUE(qvars->superclusters.empty());
    EXPECT_TRUE(tvars->superclusters.empty());
}

TEST(Supercluster, SingleCluster) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars =
            make_merged_callset(subs_at({10, 12}), {0}, {5}, {17});

    run_supercluster(qvars, make_empty_callset());

    // one cluster means one supercluster holding both variants
    EXPECT_EQ(std::vector<int>({0, 0}), qvars->superclusters);
}

TEST(Supercluster, TwoFarTwoScs) {
    GlobalsGuard guard;

    // the first cluster reaches to 15, well short of the second's leftward reach of 105
    std::shared_ptr<ctgVariants> qvars =
            make_merged_callset(subs_at({10, 110}), {0, 1}, {5, 105}, {15, 115});

    run_supercluster(qvars, make_empty_callset());

    EXPECT_EQ(std::vector<int>({0, 1}), qvars->superclusters);
}

TEST(Supercluster, QueryTruthOverlapOneSc) {
    GlobalsGuard guard;

    // the truth cluster reaches left to 12, inside the query cluster reaching right to 20
    std::shared_ptr<ctgVariants> qvars = make_merged_callset(subs_at({10}), {0}, {5}, {20});
    std::shared_ptr<ctgVariants> tvars = make_merged_callset(subs_at({14}), {0}, {12}, {25});

    run_supercluster(qvars, tvars);

    // grouping is across callsets, so both land in supercluster 0
    EXPECT_EQ(std::vector<int>({0}), qvars->superclusters);
    EXPECT_EQ(std::vector<int>({0}), tvars->superclusters);
}

TEST(Supercluster, ChainedOverlap) {
    GlobalsGuard guard;

    // query A spans [0, 20] and query C spans [40, 60], which do not touch; truth B spans
    // [15, 45] and overlaps both, so the three must chain into a single supercluster
    std::shared_ptr<ctgVariants> qvars =
            make_merged_callset(subs_at({10, 50}), {0, 1}, {0, 40}, {20, 60});
    std::shared_ptr<ctgVariants> tvars = make_merged_callset(subs_at({30}), {0}, {15}, {45});

    run_supercluster(qvars, tvars);

    // A and C are transitively joined through B, rather than splitting into two superclusters
    EXPECT_EQ(std::vector<int>({0, 0}), qvars->superclusters);
    EXPECT_EQ(std::vector<int>({0}), tvars->superclusters);
}

TEST(Supercluster, OversizedSplitInvoked) {
    GlobalsGuard guard;
    g.max_supercluster_size = 150;

    // four insertions spread over 300bp held in one cluster, so the supercluster spans 302bp
    std::shared_ptr<ctgVariants> qvars = make_merged_callset(
            {{0, 0, TYPE_INS, "", "A"}, {100, 0, TYPE_INS, "", "A"},
             {200, 0, TYPE_INS, "", "A"}, {300, 0, TYPE_INS, "", "A"}}, {0}, {-1}, {301});

    testing::internal::CaptureStderr();
    run_supercluster(qvars, make_empty_callset());
    std::string err = testing::internal::GetCapturedStderr();

    // the oversized supercluster is split at the central gap and the split is reported
    EXPECT_EQ(std::vector<int>({0, 0, 1, 1}), qvars->superclusters);
    EXPECT_NE(std::string::npos, err.find("Max supercluster size (150) exceeded (302)"));
    EXPECT_NE(std::string::npos, err.find("breaking up into 2 superclusters"));
}

TEST(Supercluster, VarCountStatGuard) {
    GlobalsGuard guard;

    // the query carries every variant; the truth callset holds none
    std::shared_ptr<ctgVariants> qvars = make_merged_callset(subs_at({10, 12}), {0}, {5}, {17});
    std::shared_ptr<ctgVariants> tvars = make_empty_callset();

    run_supercluster(qvars, tvars);

    // the per-supercluster variant-count stat is guarded by `if (vars[ci]->nc)`, so the empty
    // truth callset contributes nothing rather than reading its cluster lanes
    EXPECT_EQ(std::vector<int>({0, 0}), qvars->superclusters);
    EXPECT_EQ(0, tvars->nc);
    EXPECT_EQ(std::vector<int>({0}), tvars->clusters);

    // make_empty_callset() supplies the single boundary that the unguarded supercluster-assignment
    // loops read for every callset; a callset with a genuinely empty clusters vector would be
    // indexed out of bounds there, which is why the merge now always writes this boundary (#166)
    EXPECT_FALSE(tvars->clusters.empty());
}

TEST(Supercluster, MergedEmptyCallsetSurvives) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000},
            {make_ctgSuperclusters(make_ctgVariants("chr1", {}), make_ctgVariants("chr1", {}))});
    std::shared_ptr<ctgVariants> qhap = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_ALT}});
    set_clusters(qhap, {0, 1}, {5}, {15});
    auto qvars = make_hap_vars(qhap, make_ctgVariants("chr1", {}));
    auto tvars = make_hap_vars(make_ctgVariants("chr1", {}), make_ctgVariants("chr1", {}));

    // unlike VarCountStatGuard above, both callsets come from the merge rather than being built by
    // hand, so this composes the two functions the way the binary does: the truth callset is the
    // merge's empty-contig product, and the query has variants so the contig is not skipped
    sc_data->load_and_merge_callset_vars_across_haps(QUERY, qvars);
    sc_data->load_and_merge_callset_vars_across_haps(TRUTH, tvars);
    sc_data->supercluster();

    EXPECT_EQ(std::vector<int>({0}),
            sc_data->superclusters["chr1"]->callset_vars[QUERY]->superclusters);
    EXPECT_TRUE(sc_data->superclusters["chr1"]->callset_vars[TRUTH]->superclusters.empty());
}

TEST(Supercluster, BrksAdvanceRemainder) {
    GlobalsGuard guard;

    // an int::max left reach on the second cluster stops the main loop early while variants
    // remain, which is the only way to reach the trailing "add remaining variants" loop; a
    // sentinel-terminated callset from the merge never gets here, because its final boundary
    // already equals n
    std::shared_ptr<ctgVariants> qvars = make_ctgVariants("chr1", subs_at({10, 110}));
    set_clusters(qvars, {0, 1}, {5, INT_MAXIMUM}, {15, INT_MAXIMUM}, 2);

    run_supercluster(qvars, make_empty_callset());

    // variant 0 is assigned by the main loop, variant 1 by the remainder loop, both getting the
    // supercluster index the main loop stopped at
    EXPECT_EQ(std::vector<int>({0, 1}), qvars->superclusters);
}

/* sort_superclusters *****************************************************************************/

/**
 * @brief Builds a callset carrying per-variant supercluster assignments and a non-zero nc.
 *
 * sort_superclusters() reads only nc's non-zero-ness, so the cluster lane is nominal.
 * @param[in] vars Variants to hold, each with its supercluster field set and ascending
 * @param[in] ctg Contig name
 * @return Variant container ready to be sorted
 */
std::shared_ptr<ctgVariants> make_sorted_callset(const std::vector<var_desc> & vars,
        const std::string & ctg = "chr1") {
    std::shared_ptr<ctgVariants> sorted = make_ctgVariants(ctg, vars);
    set_clusters(sorted, {0, sorted->n}, {0}, {0});
    return sorted;
}

/**
 * @brief Describes one substitution already assigned to a supercluster.
 * @param[in] pos Reference position
 * @param[in] supercluster Supercluster index this variant belongs to
 * @param[in] gt Genotype, which decides the haplotypes the ALT length is counted on
 * @return Variant descriptor ready for make_sorted_callset()
 */
var_desc sub_in_sc(int pos, int supercluster, gt_t gt = GT_ALT_ALT) {
    var_desc var;
    var.pos = pos;
    var.rlen = 1;
    var.type = TYPE_SUB;
    var.ref = "A";
    var.alt = "C";
    var.gt = gt;
    var.supercluster = supercluster;
    return var;
}

/**
 * @brief Counts every (ctg_idx, sc_idx) pair across all buckets.
 * @param[in] groups Bucketed superclusters as returned by sort_superclusters()
 * @return Total number of superclusters placed into any bucket
 */
int total_sorted(const std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > & groups) {
    int total = 0;
    for (const auto & bucket : groups) total += int(bucket[SC_IDX].size());
    return total;
}

/*
 * These cases rely on the default thread/RAM scheduling steps, which Globals' constructor derives
 * from the default max_threads and max_ram of 8: four buckets with ceilings of 1, 2, 4 and 8 GB.
 * The estimate a supercluster is placed by is (query_len + truth_len) * 8 bytes * g.max_size * 2,
 * so with the default max_size of 1000 a single-base variant per callset needs ~32kB and a
 * 100000-base one ~1.6GB, straddling the first bucket boundary. Only a case whose assertion is
 * about the placement itself overrides a global.
 */

TEST(SortSuperclusters, EmptyQueryStillScheduled) {
    GlobalsGuard guard;

    // the truth callset holds a variant the query never called, so its supercluster still has to
    // be evaluated for that call to be counted as a false negative; the count was once keyed on
    // QUERY alone, which dropped the contig entirely (#166)
    std::shared_ptr<ctgVariants> tvars = make_sorted_callset({sub_in_sc(10, 0)});
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000}, {make_ctgSuperclusters(make_empty_callset(), tvars)});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    // one bucket per scheduling step, with the truth's lone supercluster placed in the first
    ASSERT_EQ(size_t(g.thread_nsteps), groups.size());
    EXPECT_EQ(1, total_sorted(groups));
    EXPECT_EQ(std::vector<int>({0}), groups[0][SC_IDX]);
    EXPECT_EQ(std::vector<int>({0}), groups[0][CTG_IDX]);
}

TEST(SortSuperclusters, EmptyTruthStillScheduled) {
    GlobalsGuard guard;

    // the mirror case, which read the truth's supercluster lane at index n-1 == -1 (#166)
    std::shared_ptr<ctgVariants> qvars = make_sorted_callset({sub_in_sc(10, 0)});
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000}, {make_ctgSuperclusters(qvars, make_empty_callset())});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    EXPECT_EQ(1, total_sorted(groups));
    EXPECT_EQ(std::vector<int>({0}), groups[0][SC_IDX]);
}

TEST(SortSuperclusters, BothCallsetsEmptySkipped) {
    GlobalsGuard guard;

    // a contig can reach here with no variants on either callset, when every variant on it was
    // filtered out; it holds no superclusters, so it must still contribute nothing
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000}, {make_ctgSuperclusters(make_empty_callset(), make_empty_callset())});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    ASSERT_EQ(size_t(g.thread_nsteps), groups.size());
    EXPECT_EQ(0, total_sorted(groups));
}

TEST(SortSuperclusters, NscsCount) {
    GlobalsGuard guard;

    // the query's last supercluster is 1 and the truth's is 2, so the count is max(1, 2) + 1
    std::shared_ptr<ctgVariants> qvars =
            make_sorted_callset({sub_in_sc(10, 0), sub_in_sc(20, 1)});
    std::shared_ptr<ctgVariants> tvars =
            make_sorted_callset({sub_in_sc(12, 0), sub_in_sc(30, 2)});
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000}, {make_ctgSuperclusters(qvars, tvars)});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    // supercluster 1 holds no truth variants and supercluster 2 no query variants, yet all three
    // indices are emitted
    EXPECT_EQ(3, total_sorted(groups));
    EXPECT_EQ(std::vector<int>({0, 1, 2}), groups[0][SC_IDX]);
}

TEST(SortSuperclusters, SmallLowBucket) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_sorted_callset({sub_in_sc(10, 0)});
    std::shared_ptr<ctgVariants> tvars = make_sorted_callset({sub_in_sc(12, 0)});
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000}, {make_ctgSuperclusters(qvars, tvars)});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    // ~32kB is far below the first bucket's 1GB ceiling, so nothing spills past it
    EXPECT_EQ(std::vector<int>({0}), groups[0][SC_IDX]);
    EXPECT_EQ(std::vector<int>({0}), groups[0][CTG_IDX]);
    EXPECT_EQ(1, total_sorted(groups));
}

TEST(SortSuperclusters, LargeLastBucketWarn) {
    GlobalsGuard guard;

    // a RAM ceiling below even the smallest supercluster's estimate; the placement is what this
    // case asserts, so it sets the ceiling rather than taking the default
    g.max_ram = 1e-6;
    std::shared_ptr<ctgVariants> qvars = make_sorted_callset({sub_in_sc(10, 0)});
    std::shared_ptr<ctgVariants> tvars = make_sorted_callset({sub_in_sc(12, 0)});
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000}, {make_ctgSuperclusters(qvars, tvars)});

    testing::internal::CaptureStderr();
    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);
    std::string err = testing::internal::GetCapturedStderr();

    // the supercluster is run anyway, in the last bucket, where the fewest threads are active
    EXPECT_TRUE(groups[0][SC_IDX].empty());
    EXPECT_EQ(std::vector<int>({0}), groups[g.thread_nsteps-1][SC_IDX]);
    EXPECT_NE(std::string::npos, err.find("RAM exceeded"));
    EXPECT_NE(std::string::npos, err.find("running anyways"));
}

TEST(SortSuperclusters, CtgSuperclusterPaired) {
    GlobalsGuard guard;
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1", "chr2"}, {1000, 1000},
            {make_ctgSuperclusters(make_sorted_callset({sub_in_sc(10, 0)}, "chr1"),
                                   make_sorted_callset({sub_in_sc(12, 0)}, "chr1")),
             make_ctgSuperclusters(make_sorted_callset({sub_in_sc(10, 0), sub_in_sc(20, 1)}, "chr2"),
                                   make_sorted_callset({sub_in_sc(12, 0), sub_in_sc(22, 1)}, "chr2"))});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    // the two index lanes stay parallel, so entry i names supercluster sc_idx[i] on contig
    // ctg_idx[i]; chr1 contributes one supercluster and chr2 two
    EXPECT_EQ(std::vector<int>({0, 1, 1}), groups[0][CTG_IDX]);
    EXPECT_EQ(std::vector<int>({0, 0, 1}), groups[0][SC_IDX]);
    EXPECT_EQ(groups[0][CTG_IDX].size(), groups[0][SC_IDX].size());
}

TEST(SortSuperclusters, LenLowerUpperBound) {
    GlobalsGuard guard;

    // supercluster 0 spans positions 0 to 100000, supercluster 1 holds a single variant
    std::shared_ptr<ctgVariants> qvars = make_sorted_callset(
            {sub_in_sc(0, 0), sub_in_sc(100000, 0), sub_in_sc(100010, 1)});
    std::shared_ptr<ctgVariants> tvars = make_sorted_callset(
            {sub_in_sc(50, 0), sub_in_sc(100012, 1)});
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {200000}, {make_ctgSuperclusters(qvars, tvars)});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    // the bounds select each supercluster's own variants, so the wide supercluster 0 is estimated
    // at ~1.6GB and crosses the first bucket's 1GB ceiling, while the narrow supercluster 1 stays
    // under it
    EXPECT_EQ(std::vector<int>({1}), groups[0][SC_IDX]);
    EXPECT_EQ(std::vector<int>({0}), groups[1][SC_IDX]);
}

TEST(SortSuperclusters, EmptyCallsetNcZeroGuard) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_sorted_callset({sub_in_sc(10, 0)});

    // an unclustered truth callset carrying a 100000-base alternate allele; counting it would
    // push the estimate to ~1.6GB and past the first bucket's ceiling
    var_desc big = sub_in_sc(12, 0);
    big.alt = std::string(100000, 'C');
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1", {big});
    ASSERT_EQ(0, tvars->nc);
    std::shared_ptr<superclusterData> sc_data = make_superclusterData(
            {"chr1"}, {1000}, {make_ctgSuperclusters(qvars, tvars)});

    std::vector< EnumArray<idxdim_t, std::vector<int>, IDXDIM_SLOTS> > groups = sort_superclusters(sc_data);

    // the nc == 0 guard skips that callset's length entirely, so only the query contributes
    EXPECT_EQ(std::vector<int>({0}), groups[0][SC_IDX]);
    EXPECT_TRUE(groups[1][SC_IDX].empty());
}

/* superclusterData ctor **************************************************************************/

/**
 * @brief Declares a contig on a callset without giving it any variants.
 *
 * parse_variants() creates a ctgVariants for every contig in the VCF header but only appends to
 * `contigs` once a record is seen, so a header contig with no records is reachable through
 * `variants` and absent from `contigs`. The merge indexes `variants[hap][ctg]` for every contig in
 * the union of both callsets and dereferences the result without a null check, so a contig that
 * one callset never declared at all is not a state these tests construct.
 * @param[in,out] vars Callset to declare the contig on
 * @param[in] ctg Contig name
 */
void declare_contig(std::shared_ptr<variantData> vars, const std::string & ctg) {
    vars->variants[HAP1][ctg] = make_ctgVariants(ctg, {});
    vars->variants[HAP2][ctg] = make_ctgVariants(ctg, {});
}

TEST(SuperclusterDataCtor, ContigUnionDedup) {
    GlobalsGuard guard;
    std::shared_ptr<variantData> qvd =
            make_variantData(QUERY, {"chr1", "chr2"}, {100, 200}, {{2}, {2}});
    std::shared_ptr<variantData> tvd =
            make_variantData(TRUTH, {"chr2", "chr3"}, {999, 300}, {{2}, {1}});
    declare_contig(qvd, "chr3");
    declare_contig(tvd, "chr1");

    superclusterData sc_data(qvd, tvd, nullptr);

    // query contigs come first and truth adds only what query did not already cover
    EXPECT_EQ(std::vector<std::string>({"chr1", "chr2", "chr3"}), sc_data.contigs);

    // chr2 is shared, and the query's length wins because query is scanned first
    EXPECT_EQ(std::vector<int>({100, 200, 300}), sc_data.lengths);
    EXPECT_EQ(size_t(3), sc_data.superclusters.size());
}

TEST(SuperclusterDataCtor, SamplesFilenamesOrder) {
    GlobalsGuard guard;
    std::shared_ptr<variantData> qvd = make_variantData(QUERY, {"chr1"}, {100}, {{2}});
    std::shared_ptr<variantData> tvd = make_variantData(TRUTH, {"chr1"}, {100}, {{2}});

    superclusterData sc_data(qvd, tvd, nullptr);

    // both lanes are indexed by callset, so query occupies slot QUERY and truth slot TRUTH
    ASSERT_EQ(size_t(CALLSETS), sc_data.samples.size());
    ASSERT_EQ(size_t(CALLSETS), sc_data.filenames.size());
    EXPECT_EQ("QUERY", sc_data.samples[QUERY]);
    EXPECT_EQ("TRUTH", sc_data.samples[TRUTH]);
    EXPECT_EQ("QUERY.vcf", sc_data.filenames[QUERY]);
    EXPECT_EQ("TRUTH.vcf", sc_data.filenames[TRUTH]);
}

TEST(SuperclusterDataCtor, QueryOnlyContig) {
    GlobalsGuard guard;
    std::shared_ptr<variantData> qvd = make_variantData(QUERY, {"chr1"}, {100}, {{2}});
    std::shared_ptr<variantData> tvd = make_variantData(TRUTH, {}, {}, {});
    declare_contig(tvd, "chr1");

    superclusterData sc_data(qvd, tvd, nullptr);

    // the contig survives on the query's length and ploidy, and both callsets get a container
    EXPECT_EQ(std::vector<std::string>({"chr1"}), sc_data.contigs);
    EXPECT_EQ(std::vector<int>({100}), sc_data.lengths);
    EXPECT_NE(nullptr, sc_data.superclusters["chr1"]->callset_vars[QUERY]);
    EXPECT_NE(nullptr, sc_data.superclusters["chr1"]->callset_vars[TRUTH]);
    EXPECT_EQ(0, sc_data.superclusters["chr1"]->callset_vars[TRUTH]->n);
}

TEST(SuperclusterDataCtor, TruthOnlyContig) {
    GlobalsGuard guard;
    std::shared_ptr<variantData> qvd = make_variantData(QUERY, {}, {}, {});
    std::shared_ptr<variantData> tvd = make_variantData(TRUTH, {"chr1"}, {100}, {{1}});
    declare_contig(qvd, "chr1");

    superclusterData sc_data(qvd, tvd, nullptr);

    // the truth-only pass appends the contig with the truth's length
    EXPECT_EQ(std::vector<std::string>({"chr1"}), sc_data.contigs);
    EXPECT_EQ(std::vector<int>({100}), sc_data.lengths);
    EXPECT_EQ(0, sc_data.superclusters["chr1"]->callset_vars[QUERY]->n);
}

TEST(SuperclusterDataCtor, EndToEndSmoke) {
    GlobalsGuard guard;
    std::shared_ptr<variantData> qvd = make_variantData(QUERY, {"chr1"}, {1000}, {{2}});
    std::shared_ptr<variantData> tvd = make_variantData(TRUTH, {"chr1"}, {1000}, {{2}});

    // one query variant and one nearby truth variant, whose reaches overlap, plus a distant
    // query variant that must not join them
    qvd->variants[HAP1]["chr1"] = make_ctgVariants("chr1",
            {{10, 1, TYPE_SUB, "A", "C", GT_ALT_REF},
             {500, 1, TYPE_SUB, "A", "G", GT_ALT_REF}});
    set_clusters(qvd->variants[HAP1]["chr1"], {0, 1, 2}, {5, 495}, {20, 505});
    tvd->variants[HAP2]["chr1"] = make_ctgVariants("chr1",
            {{14, 1, TYPE_SUB, "A", "T", GT_REF_ALT}});
    set_clusters(tvd->variants[HAP2]["chr1"], {0, 1}, {12}, {25});

    superclusterData sc_data(qvd, tvd, nullptr);

    // the ctor merges across haplotypes and then superclusters, so both lanes are populated
    std::shared_ptr<ctgVariants> qvars = sc_data.superclusters["chr1"]->callset_vars[QUERY];
    std::shared_ptr<ctgVariants> tvars = sc_data.superclusters["chr1"]->callset_vars[TRUTH];
    ASSERT_EQ(2, qvars->n);
    ASSERT_EQ(1, tvars->n);

    // the overlapping query and truth variants share supercluster 0, and the distant query
    // variant starts supercluster 1
    EXPECT_EQ(std::vector<int>({0, 1}), qvars->superclusters);
    EXPECT_EQ(std::vector<int>({0}), tvars->superclusters);
}

/* wf_swg_cluster *********************************************************************************/

/**
 * @brief Builds a deterministic pseudo-random reference sequence.
 *
 * A periodic reference lets an alignment slide a variant by one whole period at no cost, which
 * inflates every reach and makes cluster boundaries depend on the period rather than on the
 * penalties under test. A fixed-seed pseudo-random sequence keeps reaches tight and reproducible.
 * @param[in] length Sequence length in bases
 * @param[in] seed Seed for the linear congruential generator
 * @return Sequence of the requested length over ACGT
 */
std::string pseudo_ref(int length, uint32_t seed = 1) {
    const std::string bases = "ACGT";
    std::string seq;
    uint32_t state = seed;
    for (int i = 0; i < length; i++) {
        state = state * 1103515245u + 12345u;
        seq += bases[(state >> 16) & 3];
    }
    return seq;
}

/**
 * @brief Describes a substitution whose ALT differs from the reference base at that position.
 *
 * An ALT equal to the reference would make the variant a no-op, collapsing its alignment score
 * and therefore its reach.
 * @param[in] seq Reference sequence
 * @param[in] pos Reference position of the substitution
 * @return Variant descriptor ready for make_ctgVariants()
 */
var_desc sub_vs_ref(const std::string & seq, int pos) {
    var_desc var;
    var.pos = pos;
    var.rlen = 1;
    var.type = TYPE_SUB;
    var.ref = std::string(1, seq[pos]);
    var.alt = std::string(1, seq[pos] == 'A' ? 'T' : 'A');
    return var;
}

/**
 * @brief Builds a single-contig callset over a pseudo-random reference, ready for clustering.
 * @param[in] poss Substitution positions, ascending
 * @param[in] length Contig length
 * @param[in] hap Haplotype to place the variants on
 * @return Callset whose reference, contig length and variants are all consistent
 */
std::shared_ptr<variantData> make_cluster_input(const std::vector<int> & poss, int length = 400,
        hap_t hap = HAP1) {
    std::string seq = pseudo_ref(length);
    std::shared_ptr<variantData> vcf = make_variantData(QUERY, {"chr1"}, {length}, {{2}});
    vcf->ref = make_fasta("chr1", seq);
    std::vector<var_desc> vars;
    for (int pos : poss) vars.push_back(sub_vs_ref(seq, pos));
    vcf->variants[hap]["chr1"] = make_ctgVariants("chr1", vars);
    return vcf;
}

/**
 * @brief Sets the clustering parameters wf_swg_cluster() reads from the global `g`.
 *
 * Set explicitly rather than inherited, so these cases keep their meaning if a default moves.
 * @param[in] reach_min_gap Minimum gap bridged when merging adjacent clusters
 * @param[in] max_cluster_itrs Maximum merge-and-recompute passes
 * @param[in] max_size Maximum variant size, which sizes the reusable offsets buffer
 */
void set_cluster_params(int reach_min_gap = 10, int max_cluster_itrs = 1, int max_size = 1000) {
    g.reach_min_gap = reach_min_gap;
    g.max_cluster_itrs = max_cluster_itrs;
    g.max_size = max_size;
}

TEST(WfSwgCluster, NoVarsReturns) {
    GlobalsGuard guard;
    set_cluster_params();
    std::shared_ptr<variantData> vcf = make_cluster_input({});

    // a marker clustering that the early return must leave alone
    set_clusters(vcf->variants[HAP1]["chr1"], {7}, {8}, {9});

    wf_swg_cluster(vcf.get(), 0, HAP1, g.sub, g.open, g.extend);

    std::shared_ptr<ctgVariants> vars = vcf->variants[HAP1]["chr1"];
    EXPECT_EQ(std::vector<int>({7}), vars->clusters);
    EXPECT_EQ(std::vector<int>({8}), vars->left_reaches);
    EXPECT_EQ(std::vector<int>({9}), vars->right_reaches);
}

TEST(WfSwgCluster, SingleVariant) {
    GlobalsGuard guard;
    set_cluster_params();
    std::shared_ptr<variantData> vcf = make_cluster_input({200});

    wf_swg_cluster(vcf.get(), 0, HAP1, g.sub, g.open, g.extend);

    // one variant yields one cluster; all three lanes hold nc+1 entries, the last a sentinel
    std::shared_ptr<ctgVariants> vars = vcf->variants[HAP1]["chr1"];
    EXPECT_EQ(1, vars->nc);
    EXPECT_EQ(std::vector<int>({0, 1}), vars->clusters);
    ASSERT_EQ(size_t(2), vars->left_reaches.size());
    ASSERT_EQ(size_t(2), vars->right_reaches.size());
    EXPECT_EQ(INT_MAXIMUM, vars->left_reaches[1]);
    EXPECT_EQ(INT_MAXIMUM, vars->right_reaches[1]);

    // a lone substitution on a non-repetitive reference earns no reach beyond itself
    EXPECT_EQ(200, vars->left_reaches[0]);
    EXPECT_EQ(201, vars->right_reaches[0]);
}

TEST(WfSwgCluster, TwoAdjacentMerge) {
    GlobalsGuard guard;
    set_cluster_params();

    // two substitutions two bases apart, well inside each other's reach
    std::shared_ptr<variantData> vcf = make_cluster_input({200, 202});

    wf_swg_cluster(vcf.get(), 0, HAP1, g.sub, g.open, g.extend);

    std::shared_ptr<ctgVariants> vars = vcf->variants[HAP1]["chr1"];
    EXPECT_EQ(1, vars->nc);
    EXPECT_EQ(std::vector<int>({0, 2}), vars->clusters);
}

TEST(WfSwgCluster, TwoFarSeparate) {
    GlobalsGuard guard;
    set_cluster_params();

    // 300 bases apart, far beyond any reach a single substitution can earn
    std::shared_ptr<variantData> vcf = make_cluster_input({50, 350});

    wf_swg_cluster(vcf.get(), 0, HAP1, g.sub, g.open, g.extend);

    std::shared_ptr<ctgVariants> vars = vcf->variants[HAP1]["chr1"];
    EXPECT_EQ(2, vars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2}), vars->clusters);

    // the first cluster's reach stops short of the second cluster's
    EXPECT_LT(vars->right_reaches[0], vars->left_reaches[1]);
}

TEST(WfSwgCluster, MaxIterationsBreak) {
    GlobalsGuard guard;

    // two substitutions close enough to merge on the first pass, so the loop condition still
    // holds afterwards and only the iteration cap can stop it
    set_cluster_params(/* reach_min_gap = */ 10, /* max_cluster_itrs = */ 1);
    std::shared_ptr<variantData> one_itr = make_cluster_input({200, 202});
    wf_swg_cluster(one_itr.get(), 0, HAP1, g.sub, g.open, g.extend);

    set_cluster_params(/* reach_min_gap = */ 10, /* max_cluster_itrs = */ 4);
    std::shared_ptr<variantData> four_itrs = make_cluster_input({200, 202});
    wf_swg_cluster(four_itrs.get(), 0, HAP1, g.sub, g.open, g.extend);

    std::shared_ptr<ctgVariants> capped = one_itr->variants[HAP1]["chr1"];
    std::shared_ptr<ctgVariants> converged = four_itrs->variants[HAP1]["chr1"];

    // both agree on the clustering itself
    EXPECT_EQ(std::vector<int>({0, 2}), capped->clusters);
    EXPECT_EQ(capped->clusters, converged->clusters);

    // but the cap stops before the merged cluster's reaches are recomputed: one pass keeps the
    // minimum of the two constituent left reaches, while a second pass re-aligns the merged
    // cluster as a unit and finds it reaches one base further left
    EXPECT_EQ(200, capped->left_reaches[0]);
    EXPECT_EQ(199, converged->left_reaches[0]);
    EXPECT_EQ(203, capped->right_reaches[0]);
    EXPECT_EQ(203, converged->right_reaches[0]);
}

TEST(WfSwgCluster, ReachMinGapBoundary) {
    GlobalsGuard guard;

    // a lone substitution reaches [pos, pos+1], so cluster 0 ends at 201 and cluster 1 starts at
    // 212; the merge test is `right_reach + reach_min_gap >= left_reach`, making 11 the threshold
    set_cluster_params(/* reach_min_gap = */ 10);
    std::shared_ptr<variantData> below = make_cluster_input({200, 212});
    wf_swg_cluster(below.get(), 0, HAP1, g.sub, g.open, g.extend);

    set_cluster_params(/* reach_min_gap = */ 11);
    std::shared_ptr<variantData> at = make_cluster_input({200, 212});
    wf_swg_cluster(at.get(), 0, HAP1, g.sub, g.open, g.extend);

    // one short of the gap the clusters stay apart; exactly at it they merge
    EXPECT_EQ(2, below->variants[HAP1]["chr1"]->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2}), below->variants[HAP1]["chr1"]->clusters);
    EXPECT_EQ(1, at->variants[HAP1]["chr1"]->nc);
    EXPECT_EQ(std::vector<int>({0, 2}), at->variants[HAP1]["chr1"]->clusters);
}

TEST(WfSwgCluster, ContigEdgeClamp) {
    GlobalsGuard guard;
    set_cluster_params();

    // variants near both ends of a short contig, so the iterative-doubling windows would run off
    // each end and are clamped to the contig instead
    std::shared_ptr<variantData> vcf = make_cluster_input({5, 94}, /* length = */ 100);

    wf_swg_cluster(vcf.get(), 0, HAP1, g.sub, g.open, g.extend);

    // every real cluster's reaches stay on the contig; the trailing sentinel pair is excluded
    std::shared_ptr<ctgVariants> vars = vcf->variants[HAP1]["chr1"];
    ASSERT_EQ(2, vars->nc);
    ASSERT_EQ(size_t(3), vars->left_reaches.size());
    for (int c = 0; c < vars->nc; c++) {
        EXPECT_LE(0, vars->left_reaches[c]) << "cluster " << c;
        EXPECT_LE(vars->right_reaches[c], 100) << "cluster " << c;
    }
}

TEST(WfSwgCluster, ContigStartSub) {
    GlobalsGuard guard;
    set_cluster_params();

    // both reach windows open one base left of the first variant, which is -1 here (#167)
    std::shared_ptr<variantData> vcf = make_cluster_input({0});

    wf_swg_cluster(vcf.get(), 0, HAP1, g.sub, g.open, g.extend);

    // the span is the substituted base alone, as it is for SingleVariant at position 200
    std::shared_ptr<ctgVariants> vars = vcf->variants[HAP1]["chr1"];
    EXPECT_EQ(1, vars->nc);
    EXPECT_EQ(std::vector<int>({0, 1}), vars->clusters);
    EXPECT_EQ(std::vector<int>({0, INT_MAXIMUM}), vars->left_reaches);
    EXPECT_EQ(std::vector<int>({1, INT_MAXIMUM}), vars->right_reaches);
}

TEST(WfSwgCluster, ContigStartIndels) {
    GlobalsGuard guard;
    set_cluster_params();

    // a cluster's left reach never sits right of its own first variant, at position 0 as elsewhere
    for (const var_desc & var : std::vector<var_desc>{{0, 0, TYPE_INS, "", "GG"},
            {0, 2, TYPE_DEL, "", ""}, {0, 2, TYPE_CPX, "", "TCC"}}) {
        std::string seq = pseudo_ref(400);
        std::shared_ptr<variantData> vcf = make_variantData(QUERY, {"chr1"}, {400}, {{2}});
        vcf->ref = make_fasta("chr1", seq);
        var_desc v = var;
        if (v.rlen) v.ref = seq.substr(0, v.rlen); // a DEL/CPX ref allele must match the reference
        vcf->variants[HAP1]["chr1"] = make_ctgVariants("chr1", {v});

        wf_swg_cluster(vcf.get(), 0, HAP1, g.sub, g.open, g.extend);

        std::shared_ptr<ctgVariants> vars = vcf->variants[HAP1]["chr1"];
        EXPECT_EQ(0, vars->left_reaches[0]) << "type " << int(var.type);
        EXPECT_GT(vars->right_reaches[0], 0) << "type " << int(var.type);
    }
}

} // namespace
