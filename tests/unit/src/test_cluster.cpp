/**
 * @file test_cluster.cpp
 * @brief Unit tests for cluster.cpp: supercluster index and range arithmetic.
 */
#include <limits>
#include <memory>
#include <string>
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

    // no callset has an unprocessed variant, so the sentinel callset index is returned
    var_info info = get_next_variant_info(sc->callset_vars, {1, 1}, {1, 1});
    EXPECT_EQ(-1, info.callset_idx);
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
    EXPECT_EQ(2, split[QUERY]);
    EXPECT_EQ(0, split[TRUTH]);
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
    EXPECT_EQ(2, split[QUERY]);
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
    EXPECT_EQ(3, split[QUERY]);
    EXPECT_NE(1, split[QUERY]);
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
    EXPECT_EQ(2, split[QUERY]);
    EXPECT_NE(1, split[QUERY]);
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
    EXPECT_EQ(1, split[QUERY]);
    EXPECT_EQ(1, split[TRUTH]);
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
    std::vector< std::vector<int> > breakpoints = {{0, 0}, {3, 2}};

    // query variant 2 and truth variant 1 are both already cluster boundaries
    std::vector<int> split = split_cluster(sc->callset_vars, {2, 1}, breakpoints, 0);

    EXPECT_EQ(std::vector<int>({1, 1}), split);
    EXPECT_EQ(3, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 4, 6}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({5, 25, 45}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({25, 45, 65}), qvars->right_reaches);
    EXPECT_EQ(2, tvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2}), tvars->clusters);
    EXPECT_EQ(std::vector<int>({105, 115}), tvars->right_reaches);

    // nothing was inserted, so no later breakpoint moved
    EXPECT_EQ(std::vector<int>({0, 0}), breakpoints[0]);
    EXPECT_EQ(std::vector<int>({3, 2}), breakpoints[1]);
}

TEST(SplitCluster, MidclusterInserts) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector< std::vector<int> > breakpoints = {{0, 0}, {3, 0}};

    // query variant 3 sits inside cluster 1 (variants 2..3), so that cluster is split in two
    std::vector<int> split = split_cluster(sc->callset_vars, {3, 0}, breakpoints, 0);

    EXPECT_EQ(std::vector<int>({2, 0}), split);
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
    std::vector< std::vector<int> > breakpoints = {{0, 0}, {3, 0}};

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
    std::vector< std::vector<int> > breakpoints = {{0, 0}, {2, 0}, {3, 0}};

    split_cluster(sc->callset_vars, {3, 0}, breakpoints, 0);

    // inserting a cluster shifts every later cluster index on the split callset by one
    EXPECT_EQ(std::vector<int>({0, 0}), breakpoints[0]);
    EXPECT_EQ(std::vector<int>({3, 0}), breakpoints[1]);
    EXPECT_EQ(std::vector<int>({4, 0}), breakpoints[2]);

    // only breakpoints strictly after breakpoint_idx move
    std::shared_ptr<ctgVariants> qvars2 = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc2 = make_ctgSuperclusters(qvars2, make_empty_callset());
    std::vector< std::vector<int> > breakpoints2 = {{0, 0}, {2, 0}, {3, 0}};
    split_cluster(sc2->callset_vars, {3, 0}, breakpoints2, 1);
    EXPECT_EQ(std::vector<int>({2, 0}), breakpoints2[1]);
    EXPECT_EQ(std::vector<int>({4, 0}), breakpoints2[2]);
}

TEST(SplitCluster, BothCallsetsIndependent) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgVariants> tvars = make_ctgVariants("chr1",
            {{100, 1, TYPE_SUB, "A", "C"}, {110, 1, TYPE_SUB, "A", "C"}});
    set_clusters(tvars, {0, 2}, {95}, {115});
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, tvars);
    std::vector< std::vector<int> > breakpoints = {{0, 0}, {3, 1}};

    // both callsets split mid-cluster, at different cluster indices
    std::vector<int> split = split_cluster(sc->callset_vars, {3, 1}, breakpoints, 0);

    EXPECT_EQ(std::vector<int>({2, 1}), split);
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 3, 4, 6}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({5, 25, 40, 45}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({25, 40, 45, 65}), qvars->right_reaches);
    EXPECT_EQ(2, tvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2}), tvars->clusters);
    EXPECT_EQ(std::vector<int>({95, 110}), tvars->left_reaches);
    EXPECT_EQ(std::vector<int>({110, 115}), tvars->right_reaches);

    // each callset's later breakpoints are incremented independently
    EXPECT_EQ(std::vector<int>({4, 2}), breakpoints[1]);
}

TEST(SplitCluster, ClustIdxMinusOneSafety) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = make_split_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector< std::vector<int> > breakpoints = {{0, 0}, {3, 0}};

    // variant 0 equals clusters[0], so the equality branch is taken and right_reaches[clust_idx-1]
    // is never indexed at -1; callers only pass indices at or after the first boundary
    std::vector<int> split = split_cluster(sc->callset_vars, {0, 0}, breakpoints, 0);

    EXPECT_EQ(std::vector<int>({0, 0}), split);
    EXPECT_EQ(3, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 4, 6}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({5, 25, 45}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({25, 45, 65}), qvars->right_reaches);
    EXPECT_EQ(std::vector<int>({3, 0}), breakpoints[1]);
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
    std::vector<int> end_indices = {4, 0};

    // the supercluster spans 302bp, well within the limit
    std::vector< std::vector<int> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    ASSERT_EQ(size_t(2), breakpoints.size());
    EXPECT_EQ(std::vector<int>({0, 0}), breakpoints[0]);
    EXPECT_EQ(std::vector<int>({4, 0}), breakpoints[1]);
    EXPECT_EQ(std::vector<int>({4, 0}), end_indices);
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2, 3, 4}), qvars->clusters);
}

TEST(SplitLargeSupercluster, OneSplit) {
    GlobalsGuard guard;
    g.max_supercluster_size = 150;
    std::shared_ptr<ctgVariants> qvars = make_spread_query(0);
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector<int> end_indices = {4, 0};

    // 302bp exceeds the limit, and the central gap splits it into two 102bp pieces
    std::vector< std::vector<int> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    ASSERT_EQ(size_t(3), breakpoints.size());
    EXPECT_EQ(std::vector<int>({0, 0}), breakpoints[0]);
    EXPECT_EQ(std::vector<int>({2, 0}), breakpoints[1]);
    EXPECT_EQ(std::vector<int>({4, 0}), breakpoints[2]);

    // the split fell on an existing cluster boundary, so no cluster was inserted
    EXPECT_EQ(4, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 1, 2, 3, 4}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({4, 0}), end_indices);
}

TEST(SplitLargeSupercluster, MultipleSplits) {
    GlobalsGuard guard;
    g.max_supercluster_size = 80;
    std::shared_ptr<ctgVariants> qvars = make_spread_query(0);
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector<int> end_indices = {4, 0};

    // each 102bp half still exceeds 80bp, so a second round splits both of them
    std::vector< std::vector<int> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    ASSERT_EQ(size_t(5), breakpoints.size());
    for (size_t i = 0; i < breakpoints.size(); i++) {
        EXPECT_EQ(std::vector<int>({int(i), 0}), breakpoints[i]) << "breakpoint " << i;
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
    std::vector<int> end_indices = {4, 0};

    testing::internal::CaptureStderr();
    std::vector< std::vector<int> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);
    std::string err = testing::internal::GetCapturedStderr();

    // the oversized supercluster is retained unchanged
    ASSERT_EQ(size_t(2), breakpoints.size());
    EXPECT_EQ(std::vector<int>({0, 0}), breakpoints[0]);
    EXPECT_EQ(std::vector<int>({4, 0}), breakpoints[1]);
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
    std::vector<int> end_indices = {1, 0};

    std::vector< std::vector<int> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    // the in/out end indices are overwritten with the final breakpoint, which the insertion moved
    EXPECT_EQ(std::vector<int>({2, 0}), end_indices);
    EXPECT_EQ(breakpoints[breakpoints.size()-1], end_indices);
    EXPECT_NE(std::vector<int>({1, 0}), end_indices);
}

TEST(SplitLargeSupercluster, BreakpointShiftAfterInsert) {
    GlobalsGuard guard;
    g.max_supercluster_size = 150;
    std::shared_ptr<ctgVariants> qvars = make_one_cluster_query();
    std::shared_ptr<ctgSuperclusters> sc = make_ctgSuperclusters(qvars, make_empty_callset());
    std::vector<int> end_indices = {1, 0};

    std::vector< std::vector<int> > breakpoints =
            split_large_supercluster(sc->callset_vars, {0, 0}, end_indices);

    // splitting at variant 2 inserted a cluster, so the trailing breakpoint shifted 1 -> 2
    ASSERT_EQ(size_t(3), breakpoints.size());
    EXPECT_EQ(std::vector<int>({0, 0}), breakpoints[0]);
    EXPECT_EQ(std::vector<int>({1, 0}), breakpoints[1]);
    EXPECT_EQ(std::vector<int>({2, 0}), breakpoints[2]);

    // the full post-state of the mutated cluster lanes
    EXPECT_EQ(2, qvars->nc);
    EXPECT_EQ(std::vector<int>({0, 2, 4}), qvars->clusters);
    EXPECT_EQ(std::vector<int>({-1, 200}), qvars->left_reaches);
    EXPECT_EQ(std::vector<int>({200, 301}), qvars->right_reaches);
}

} // namespace
