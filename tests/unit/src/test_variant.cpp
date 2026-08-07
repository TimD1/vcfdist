/**
 * @file test_variant.cpp
 * @brief Unit tests for variant.cpp: genotype, allele-count, and variant-type logic.
 */
#include <memory>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/defs.h"
#include "../../../src/globals.h"
#include "../../../src/variant.h"
#include "test_helpers.h"

namespace {

/* ctgVariants constructor ************************************************************************/

TEST(CtgVariantsCtor, SetsCtg) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    EXPECT_EQ("chr20", vars.ctg);
}

TEST(CtgVariantsCtor, ZeroesN) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    EXPECT_EQ(0, vars.n);
    EXPECT_EQ(0, vars.nc);

    // the parsed-data and phasing lanes are sized by add_var(), so they start empty
    EXPECT_TRUE(vars.poss.empty());
    EXPECT_TRUE(vars.rlens.empty());
    EXPECT_TRUE(vars.types.empty());
    EXPECT_TRUE(vars.locs.empty());
    EXPECT_TRUE(vars.refs.empty());
    EXPECT_TRUE(vars.alts.empty());
    EXPECT_TRUE(vars.orig_gts.empty());
    EXPECT_TRUE(vars.gt_quals.empty());
    EXPECT_TRUE(vars.var_quals.empty());
    EXPECT_TRUE(vars.phase_sets.empty());
    EXPECT_TRUE(vars.rec_idxs.empty());
    EXPECT_TRUE(vars.alt_idxs.empty());
    EXPECT_TRUE(vars.ploidies.empty());
    EXPECT_TRUE(vars.superclusters.empty());
    EXPECT_TRUE(vars.matched_gts.empty());
    EXPECT_TRUE(vars.phases.empty());
    EXPECT_TRUE(vars.pb_phases.empty());
    EXPECT_TRUE(vars.ac_errtype.empty());
    EXPECT_TRUE(vars.clusters.empty());
    EXPECT_TRUE(vars.left_reaches.empty());
    EXPECT_TRUE(vars.right_reaches.empty());
}

TEST(CtgVariantsCtor, AllocatesTwoPhaseLanes) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");

    // all six per-haplotype lanes are indexed by HAP1/HAP2, so the outer vectors are sized PHASES
    ASSERT_EQ(HAP_SLOTS, vars.errtypes.size());
    ASSERT_EQ(HAP_SLOTS, vars.sync_group.size());
    ASSERT_EQ(HAP_SLOTS, vars.callq.size());
    ASSERT_EQ(HAP_SLOTS, vars.credit.size());
    ASSERT_EQ(HAP_SLOTS, vars.ref_ed.size());
    ASSERT_EQ(HAP_SLOTS, vars.query_ed.size());

    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
        EXPECT_TRUE(vars.errtypes[hap].empty()) << "hap " << idx(hap);
        EXPECT_TRUE(vars.sync_group[hap].empty()) << "hap " << idx(hap);
        EXPECT_TRUE(vars.callq[hap].empty()) << "hap " << idx(hap);
        EXPECT_TRUE(vars.credit[hap].empty()) << "hap " << idx(hap);
        EXPECT_TRUE(vars.ref_ed[hap].empty()) << "hap " << idx(hap);
        EXPECT_TRUE(vars.query_ed[hap].empty()) << "hap " << idx(hap);
    }
}

TEST(CtgVariantsCtor, EmptyCtg) {
    GlobalsGuard guard;
    ctgVariants vars("");
    EXPECT_EQ("", vars.ctg);
    EXPECT_EQ(0, vars.n);
    EXPECT_EQ(size_t(PHASES), vars.errtypes.size());
}

/* add_var ****************************************************************************************/

// Asserts every per-variant field that exists today, so it must gain assertions as new per-variant
// vectors land, or "every field" quietly stops being every field: strata_bits (#47), is_phased
// (#46). rec_idxs/alt_idxs/ploidies (#48) are covered below. A new var_fields member with a
// default reaches its vector only if asserted here; the compiler will not require it at any call
// site.
TEST(AddVar, AllFields) {
    GlobalsGuard guard;
    g.max_qual = 100;
    ctgVariants vars("chr20");

    // every per-haplotype field differs between haplotypes, so a hap[HAP1]/hap[HAP2] mix-up in
    // add_var()'s body cannot pass
    vars.add_var(var_fields{.pos = 500, .rlen = 2, .type = TYPE_CPX, .loc = BED_OUTSIDE, .ref = "AC",
            .alt = "GT", .orig_gt = GT_ALT1_ALT1, .gt_qual = 21, .var_qual = 22, .phase_set = 33,
            .rec_idx = 12, .alt_idx = 3, .ploidy = 2, .supercluster = 7, .matched_gt = GT_ALT1_REF,
            .hap = {{{{.errtype = ERRTYPE_FN, .sync_group = 4, .callq = 6.5, .ref_ed = 8,
                       .query_ed = 10, .credit = 0.4},
                      {.errtype = ERRTYPE_TP, .sync_group = 5, .callq = 7.5, .ref_ed = 9,
                       .query_ed = 11, .credit = 0.6}}}}});

    ASSERT_EQ(1, vars.n);
    EXPECT_EQ(500, vars.poss[0]);
    EXPECT_EQ(2, vars.rlens[0]);
    EXPECT_EQ(TYPE_CPX, vars.types[0]);
    EXPECT_EQ(BED_OUTSIDE, vars.locs[0]);
    EXPECT_EQ("AC", vars.refs[0]);
    EXPECT_EQ("GT", vars.alts[0]);
    EXPECT_EQ(GT_ALT1_ALT1, vars.orig_gts[0]);
    EXPECT_FLOAT_EQ(21, vars.gt_quals[0]);
    EXPECT_FLOAT_EQ(22, vars.var_quals[0]);
    EXPECT_EQ(33, vars.phase_sets[0]);
    EXPECT_EQ(12, vars.rec_idxs[0]);
    EXPECT_EQ(3, vars.alt_idxs[0]);
    EXPECT_EQ(2, vars.ploidies[0]);
    EXPECT_EQ(7, vars.superclusters[0]);
    EXPECT_EQ(GT_ALT1_REF, vars.matched_gts[0]);
    EXPECT_EQ(ERRTYPE_FN, vars.errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_TP, vars.errtypes[HAP2][0]);
    EXPECT_EQ(4, vars.sync_group[HAP1][0]);
    EXPECT_EQ(5, vars.sync_group[HAP2][0]);
    EXPECT_FLOAT_EQ(6.5, vars.callq[HAP1][0]);
    EXPECT_FLOAT_EQ(7.5, vars.callq[HAP2][0]);
    EXPECT_EQ(8, vars.ref_ed[HAP1][0]);
    EXPECT_EQ(9, vars.ref_ed[HAP2][0]);
    EXPECT_EQ(10, vars.query_ed[HAP1][0]);
    EXPECT_EQ(11, vars.query_ed[HAP2][0]);
    EXPECT_FLOAT_EQ(0.4, vars.credit[HAP1][0]);
    EXPECT_FLOAT_EQ(0.6, vars.credit[HAP2][0]);
}

TEST(AddVar, QualCapped) {
    GlobalsGuard guard;
    g.max_qual = 60;
    ctgVariants vars("chr20");

    // the two quals differ, so clamping the wrong field cannot pass
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 99, .var_qual = 80, .phase_set = 0});

    EXPECT_FLOAT_EQ(60, vars.gt_quals[0]);
    EXPECT_FLOAT_EQ(60, vars.var_quals[0]);
}

TEST(AddVar, GtQualCappedIndependentlyOfVarQual) {
    GlobalsGuard guard;
    g.max_qual = 60;
    ctgVariants vars("chr20");

    // only gt_qual exceeds the cap, so its clamp cannot be riding on var_qual's
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 99, .var_qual = 30, .phase_set = 0});

    EXPECT_FLOAT_EQ(60, vars.gt_quals[0]);
    EXPECT_FLOAT_EQ(30, vars.var_quals[0]);
}

TEST(AddVar, QualBelowCap) {
    GlobalsGuard guard;
    g.max_qual = 60;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});
    EXPECT_FLOAT_EQ(30, vars.gt_quals[0]);
    EXPECT_FLOAT_EQ(30, vars.var_quals[0]);
}

TEST(AddVar, QualNegative) {
    GlobalsGuard guard;
    g.max_qual = 60;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = -5, .var_qual = -5, .phase_set = 0});

    // std::min() only caps from above, so a negative quality is stored unchanged
    EXPECT_FLOAT_EQ(-5, vars.gt_quals[0]);
    EXPECT_FLOAT_EQ(-5, vars.var_quals[0]);
}

TEST(AddVar, PushesPhaseDefaults) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});

    ASSERT_EQ(size_t(1), vars.phases.size());
    ASSERT_EQ(size_t(1), vars.pb_phases.size());
    ASSERT_EQ(size_t(1), vars.ac_errtype.size());
    EXPECT_EQ(PHASE_NONE, vars.phases[0]);
    EXPECT_EQ(PHASE_NONE, vars.pb_phases[0]);
    EXPECT_EQ(AC_UNKNOWN, vars.ac_errtype[0]);
}

TEST(AddVar, HeaderDefaults) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});

    // the rec_idx/alt_idx/ploidy defaults are asserted by ProvenanceVectorDefaults.UnknownSentinels
    EXPECT_EQ(-1, vars.superclusters[0]);
    EXPECT_EQ(GT_REF_REF, vars.matched_gts[0]);
    EXPECT_EQ(ERRTYPE_UN, vars.errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_UN, vars.errtypes[HAP2][0]);
    EXPECT_EQ(0, vars.sync_group[HAP1][0]);
    EXPECT_EQ(0, vars.sync_group[HAP2][0]);
    EXPECT_FLOAT_EQ(0, vars.callq[HAP1][0]);
    EXPECT_FLOAT_EQ(0, vars.callq[HAP2][0]);
    EXPECT_EQ(0, vars.ref_ed[HAP1][0]);
    EXPECT_EQ(0, vars.ref_ed[HAP2][0]);
    EXPECT_EQ(0, vars.query_ed[HAP1][0]);
    EXPECT_EQ(0, vars.query_ed[HAP2][0]);
    EXPECT_FLOAT_EQ(0, vars.credit[HAP1][0]);
    EXPECT_FLOAT_EQ(0, vars.credit[HAP2][0]);
}

TEST(AddVar, AppendsNotOverwrites) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 50, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});
    vars.add_var(var_fields{.pos = 300, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "G",
            .alt = "T", .orig_gt = GT_ALT1_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});

    ASSERT_EQ(2, vars.n);
    EXPECT_EQ(50, vars.poss[0]);
    EXPECT_EQ("A", vars.refs[0]);
    EXPECT_EQ(GT_REF_ALT1, vars.orig_gts[0]);
    EXPECT_EQ(300, vars.poss[1]);
    EXPECT_EQ("G", vars.refs[1]);
    EXPECT_EQ(GT_ALT1_ALT1, vars.orig_gts[1]);
}

TEST(AddVar, LaneLengthsTrackN) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    for (int i = 0; i < 3; i++) {
        vars.add_var(var_fields{.pos = 100*i, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
                .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});
    }

    ASSERT_EQ(3, vars.n);
    size_t n = size_t(vars.n);
    EXPECT_EQ(n, vars.poss.size());
    EXPECT_EQ(n, vars.rlens.size());
    EXPECT_EQ(n, vars.types.size());
    EXPECT_EQ(n, vars.locs.size());
    EXPECT_EQ(n, vars.refs.size());
    EXPECT_EQ(n, vars.alts.size());
    EXPECT_EQ(n, vars.orig_gts.size());
    EXPECT_EQ(n, vars.gt_quals.size());
    EXPECT_EQ(n, vars.var_quals.size());
    EXPECT_EQ(n, vars.phase_sets.size());
    EXPECT_EQ(n, vars.rec_idxs.size());
    EXPECT_EQ(n, vars.alt_idxs.size());
    EXPECT_EQ(n, vars.ploidies.size());
    EXPECT_EQ(n, vars.superclusters.size());
    EXPECT_EQ(n, vars.matched_gts.size());
    EXPECT_EQ(n, vars.phases.size());
    EXPECT_EQ(n, vars.pb_phases.size());
    EXPECT_EQ(n, vars.ac_errtype.size());
    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
        EXPECT_EQ(n, vars.errtypes[hap].size()) << "hap " << idx(hap);
        EXPECT_EQ(n, vars.sync_group[hap].size()) << "hap " << idx(hap);
        EXPECT_EQ(n, vars.callq[hap].size()) << "hap " << idx(hap);
        EXPECT_EQ(n, vars.ref_ed[hap].size()) << "hap " << idx(hap);
        EXPECT_EQ(n, vars.query_ed[hap].size()) << "hap " << idx(hap);
        EXPECT_EQ(n, vars.credit[hap].size()) << "hap " << idx(hap);
    }
}

TEST(AddVar, InsRlenZero) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 0, .type = TYPE_INS, .loc = BED_INSIDE, .ref = "",
            .alt = "ACGT", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});

    ASSERT_EQ(1, vars.n);
    EXPECT_EQ(0, vars.rlens[0]);
    EXPECT_EQ(TYPE_INS, vars.types[0]);
    EXPECT_EQ("", vars.refs[0]);
    EXPECT_EQ("ACGT", vars.alts[0]);
}

TEST(AddVar, DelEmptyAlt) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 4, .type = TYPE_DEL, .loc = BED_INSIDE, .ref = "ACGT",
            .alt = "", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});

    ASSERT_EQ(1, vars.n);
    EXPECT_EQ(4, vars.rlens[0]);
    EXPECT_EQ(TYPE_DEL, vars.types[0]);
    EXPECT_EQ("ACGT", vars.refs[0]);
    EXPECT_EQ("", vars.alts[0]);
}

// Pins every optional field's default, so a default that drifts fails here rather than silently
// changing what an omitted field means at the 20-odd call sites that rely on it.
TEST(AddVar, OptionalFieldDefaults) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});

    ASSERT_EQ(1, vars.n);
    EXPECT_EQ(-1, vars.rec_idxs[0]);
    EXPECT_EQ(-1, vars.alt_idxs[0]);
    EXPECT_EQ(0, vars.ploidies[0]);
    EXPECT_EQ(-1, vars.superclusters[0]);
    EXPECT_EQ(GT_REF_REF, vars.matched_gts[0]);
    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
        EXPECT_EQ(ERRTYPE_UN, vars.errtypes[hap][0]) << "hap " << idx(hap);
        EXPECT_EQ(0, vars.sync_group[hap][0]) << "hap " << idx(hap);
        EXPECT_FLOAT_EQ(0, vars.callq[hap][0]) << "hap " << idx(hap);
        EXPECT_EQ(0, vars.ref_ed[hap][0]) << "hap " << idx(hap);
        EXPECT_EQ(0, vars.query_ed[hap][0]) << "hap " << idx(hap);
        EXPECT_FLOAT_EQ(0, vars.credit[hap][0]) << "hap " << idx(hap);
    }
}

/* get_var ****************************************************************************************/

// get_var() exists so callers copying a variant between containers do not hand-copy parallel
// vectors, so the round trip must carry every field: any omission reintroduces the miswiring the
// accessor was added to prevent. Per-haplotype values differ so a transposition cannot pass.
TEST(GetVar, RoundTripsEveryField) {
    GlobalsGuard guard;
    g.max_qual = 100;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 500, .rlen = 2, .type = TYPE_CPX, .loc = BED_OUTSIDE, .ref = "AC",
            .alt = "GT", .orig_gt = GT_ALT1_ALT1, .gt_qual = 21, .var_qual = 22, .phase_set = 33,
            .rec_idx = 12, .alt_idx = 3, .ploidy = 2, .supercluster = 7, .matched_gt = GT_ALT1_REF,
            .hap = {{{{.errtype = ERRTYPE_FN, .sync_group = 4, .callq = 6.5, .ref_ed = 8,
                       .query_ed = 10, .credit = 0.4},
                      {.errtype = ERRTYPE_TP, .sync_group = 5, .callq = 7.5, .ref_ed = 9,
                       .query_ed = 11, .credit = 0.6}}}}});

    var_fields var = vars.get_var(0);
    EXPECT_EQ(500, var.pos);
    EXPECT_EQ(2, var.rlen);
    EXPECT_EQ(TYPE_CPX, var.type);
    EXPECT_EQ(BED_OUTSIDE, var.loc);
    EXPECT_EQ("AC", var.ref);
    EXPECT_EQ("GT", var.alt);
    EXPECT_EQ(GT_ALT1_ALT1, var.orig_gt);
    EXPECT_FLOAT_EQ(21, var.gt_qual);
    EXPECT_FLOAT_EQ(22, var.var_qual);
    EXPECT_EQ(33, var.phase_set);
    EXPECT_EQ(12, var.rec_idx);
    EXPECT_EQ(3, var.alt_idx);
    EXPECT_EQ(2, var.ploidy);
    EXPECT_EQ(7, var.supercluster);
    EXPECT_EQ(GT_ALT1_REF, var.matched_gt);
    EXPECT_EQ(ERRTYPE_FN, var.hap[HAP1].errtype);
    EXPECT_EQ(ERRTYPE_TP, var.hap[HAP2].errtype);
    EXPECT_EQ(4, var.hap[HAP1].sync_group);
    EXPECT_EQ(5, var.hap[HAP2].sync_group);
    EXPECT_FLOAT_EQ(6.5, var.hap[HAP1].callq);
    EXPECT_FLOAT_EQ(7.5, var.hap[HAP2].callq);
    EXPECT_EQ(8, var.hap[HAP1].ref_ed);
    EXPECT_EQ(9, var.hap[HAP2].ref_ed);
    EXPECT_EQ(10, var.hap[HAP1].query_ed);
    EXPECT_EQ(11, var.hap[HAP2].query_ed);
    EXPECT_FLOAT_EQ(0.4, var.hap[HAP1].credit);
    EXPECT_FLOAT_EQ(0.6, var.hap[HAP2].credit);
}

// The cluster.cpp merge loop reads with get_var() and appends with add_var(), so a field dropped
// by either one would be lost in the merged callset without any diagnostic.
TEST(GetVar, FeedsAddVarWithoutLoss) {
    GlobalsGuard guard;
    g.max_qual = 100;
    ctgVariants src("chr20");
    src.add_var(var_fields{.pos = 500, .rlen = 2, .type = TYPE_DEL, .loc = BED_BORDER, .ref = "AC",
            .alt = "", .orig_gt = GT_REF_ALT1, .gt_qual = 21, .var_qual = 22, .phase_set = 33,
            .rec_idx = 12, .alt_idx = 3, .ploidy = 2, .supercluster = 7, .matched_gt = GT_ALT1_REF,
            .hap = {{{{.errtype = ERRTYPE_FN, .sync_group = 4, .callq = 6.5, .ref_ed = 8,
                       .query_ed = 10, .credit = 0.4},
                      {.errtype = ERRTYPE_TP, .sync_group = 5, .callq = 7.5, .ref_ed = 9,
                       .query_ed = 11, .credit = 0.6}}}}});

    ctgVariants dst("chr20");
    dst.add_var(src.get_var(0));

    ASSERT_EQ(1, dst.n);
    EXPECT_EQ(src.poss[0], dst.poss[0]);
    EXPECT_EQ(src.rlens[0], dst.rlens[0]);
    EXPECT_EQ(src.types[0], dst.types[0]);
    EXPECT_EQ(src.locs[0], dst.locs[0]);
    EXPECT_EQ(src.refs[0], dst.refs[0]);
    EXPECT_EQ(src.alts[0], dst.alts[0]);
    EXPECT_EQ(src.orig_gts[0], dst.orig_gts[0]);
    EXPECT_FLOAT_EQ(src.gt_quals[0], dst.gt_quals[0]);
    EXPECT_FLOAT_EQ(src.var_quals[0], dst.var_quals[0]);
    EXPECT_EQ(src.phase_sets[0], dst.phase_sets[0]);
    EXPECT_EQ(src.rec_idxs[0], dst.rec_idxs[0]);
    EXPECT_EQ(src.alt_idxs[0], dst.alt_idxs[0]);
    EXPECT_EQ(src.ploidies[0], dst.ploidies[0]);
    EXPECT_EQ(src.superclusters[0], dst.superclusters[0]);
    EXPECT_EQ(src.matched_gts[0], dst.matched_gts[0]);
    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
        EXPECT_EQ(src.errtypes[hap][0], dst.errtypes[hap][0]) << "hap " << idx(hap);
        EXPECT_EQ(src.sync_group[hap][0], dst.sync_group[hap][0]) << "hap " << idx(hap);
        EXPECT_FLOAT_EQ(src.callq[hap][0], dst.callq[hap][0]) << "hap " << idx(hap);
        EXPECT_EQ(src.ref_ed[hap][0], dst.ref_ed[hap][0]) << "hap " << idx(hap);
        EXPECT_EQ(src.query_ed[hap][0], dst.query_ed[hap][0]) << "hap " << idx(hap);
        EXPECT_FLOAT_EQ(src.credit[hap][0], dst.credit[hap][0]) << "hap " << idx(hap);
    }
}

TEST(GetVar, RejectsOutOfRangeIndex) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(var_fields{.pos = 100, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "C", .orig_gt = GT_REF_ALT1, .gt_qual = 30, .var_qual = 30, .phase_set = 0});

    EXPECT_EXIT(vars.get_var(-1), testing::ExitedWithCode(1), "out of range");
    EXPECT_EXIT(vars.get_var(1), testing::ExitedWithCode(1), "out of range");
}

/* provenance and ploidy vectors ******************************************************************/

// Record ordinals within FIXTURE_RECORDS, pinned by the rec_idxs assertions below.
const int REC_HOM_SNP    = 0;
const int REC_HET_SNP    = 1;
const int REC_MULTIALLIC = 2;
const int REC_CPX        = 3;
const int REC_HAPLOID    = 4;

// One record per provenance case, in file order; ctgX carries the haploid call so that its
// ploidy of 1 does not conflict with the diploid ploidy recorded for ctg1.
const std::vector<std::string> FIXTURE_RECORDS = {
    "ctg1\t11\t.\tA\tG\t30\tPASS\t.\tGT\t1|1",       // hom SNP, both haps, ALT 1
    "ctg1\t21\t.\tC\tT\t30\tPASS\t.\tGT\t0|1",       // het SNP, HAP2 only, ALT 1
    "ctg1\t31\t.\tA\tG,T\t30\tPASS\t.\tGT\t1|2",     // multiallelic, ALT 1 HAP1, ALT 2 HAP2
    "ctg1\t41\t.\tATTT\tG,GG\t30\tPASS\t.\tGT\t0|2", // CPX from ALT 2, split into INS + DEL
    "ctgX\t11\t.\tA\tC\t30\tPASS\t.\tGT\t1",         // haploid SNP, HAP1 only, ALT 1
};

/**
 * @class ProvenanceVectors
 * @brief Parses the shared VCF fixture once per test, with global state saved and restored.
 */
class ProvenanceVectors : public ::testing::Test {
protected:
    /** @brief Writes and parses the shared fixture into the QUERY callset. */
    void SetUp() override {
        vcf_opts opts = make_vcf_opts(QUERY, {"ctg1", "ctgX"}, 200);
        std::string vcf_fn = write_tmp_vcf(dir, FIXTURE_RECORDS, opts);

        vcf_data = std::shared_ptr<variantData>(new variantData());
        // parse_variants only stores the reference pointer, so no FASTA fixture is needed
        parse_variants(vcf_fn, vcf_data, std::shared_ptr<fastaData>(nullptr), QUERY);

        hap1 = vcf_data->variants[HAP1]["ctg1"];
        hap2 = vcf_data->variants[HAP2]["ctg1"];
        hapx = vcf_data->variants[HAP1]["ctgX"];
    }

    GlobalsGuard guard;                    ///< Saves global state and silences parser output
    TempDir dir;                           ///< Holds the fixture VCF for the test's lifetime
    std::shared_ptr<variantData> vcf_data; ///< Parsed fixture
    std::shared_ptr<ctgVariants> hap1;     ///< ctg1 haplotype 1 variants
    std::shared_ptr<ctgVariants> hap2;     ///< ctg1 haplotype 2 variants
    std::shared_ptr<ctgVariants> hapx;     ///< ctgX haplotype 1 variants
};

// The three new vectors must stay sized n, like every other parsed-data vector.
TEST_F(ProvenanceVectors, VectorsSizedN) {
    for (const auto & vars : {hap1, hap2, hapx}) {
        EXPECT_EQ(size_t(vars->n), vars->rec_idxs.size());
        EXPECT_EQ(size_t(vars->n), vars->alt_idxs.size());
        EXPECT_EQ(size_t(vars->n), vars->ploidies.size());
    }
    EXPECT_EQ(2, hap1->n); // hom SNP, multiallelic ALT 1
    EXPECT_EQ(5, hap2->n); // hom SNP, het SNP, multiallelic ALT 2, CPX INS, CPX DEL
    EXPECT_EQ(1, hapx->n); // haploid SNP
    EXPECT_EQ(0, vcf_data->variants[HAP2]["ctgX"]->n);
}

// A 1|1 record yields one variant per haplotype, both from record 0's first ALT.
TEST_F(ProvenanceVectors, HomozygousSnp) {
    ASSERT_LE(1, hap1->n);
    EXPECT_EQ(10, hap1->poss[0]);
    EXPECT_EQ(REC_HOM_SNP, hap1->rec_idxs[0]);
    EXPECT_EQ(1, hap1->alt_idxs[0]);
    EXPECT_EQ(2, hap1->ploidies[0]);

    ASSERT_LE(1, hap2->n);
    EXPECT_EQ(10, hap2->poss[0]);
    EXPECT_EQ(REC_HOM_SNP, hap2->rec_idxs[0]);
    EXPECT_EQ(1, hap2->alt_idxs[0]);
    EXPECT_EQ(2, hap2->ploidies[0]);
}

// A 0|1 record yields a single HAP2 variant, carrying record 1's ordinal.
TEST_F(ProvenanceVectors, HeterozygousSnp) {
    ASSERT_LE(2, hap2->n);
    EXPECT_EQ(20, hap2->poss[1]);
    EXPECT_EQ(REC_HET_SNP, hap2->rec_idxs[1]);
    EXPECT_EQ(1, hap2->alt_idxs[1]);
    EXPECT_EQ(2, hap2->ploidies[1]);
}

// A 1|2 record splits across haplotypes, each half keeping its own original ALT ordinal.
TEST_F(ProvenanceVectors, MultiallelicRecord) {
    ASSERT_LE(2, hap1->n);
    EXPECT_EQ(30, hap1->poss[1]);
    EXPECT_EQ("G", hap1->alts[1]);
    EXPECT_EQ(REC_MULTIALLIC, hap1->rec_idxs[1]);
    EXPECT_EQ(1, hap1->alt_idxs[1]);
    EXPECT_EQ(2, hap1->ploidies[1]);

    ASSERT_LE(3, hap2->n);
    EXPECT_EQ(30, hap2->poss[2]);
    EXPECT_EQ("T", hap2->alts[2]);
    EXPECT_EQ(REC_MULTIALLIC, hap2->rec_idxs[2]);
    EXPECT_EQ(2, hap2->alt_idxs[2]);
    EXPECT_EQ(2, hap2->ploidies[2]);
}

// A single-allele record records ploidy 1, distinct from the diploid records above.
TEST_F(ProvenanceVectors, HaploidRecord) {
    ASSERT_LE(1, hapx->n);
    EXPECT_EQ(10, hapx->poss[0]);
    EXPECT_EQ(REC_HAPLOID, hapx->rec_idxs[0]);
    EXPECT_EQ(1, hapx->alt_idxs[0]);
    EXPECT_EQ(1, hapx->ploidies[0]);
}

// The INS and DEL halves of a CPX allele derive from one original allele, so they must agree on
// rec_idx and alt_idx (here ALT 2, not ALT 1) as well as position.
TEST_F(ProvenanceVectors, ComplexVariantHalvesShareAltIdx) {
    ASSERT_LE(5, hap2->n);
    const int ins = 3;
    const int del = 4;

    EXPECT_EQ(TYPE_INS, hap2->types[ins]);
    EXPECT_EQ(TYPE_DEL, hap2->types[del]);
    EXPECT_EQ(40, hap2->poss[ins]);
    EXPECT_EQ(40, hap2->poss[del]);
    EXPECT_EQ("GG", hap2->alts[ins]);
    EXPECT_EQ("ATTT", hap2->refs[del]);

    EXPECT_EQ(2, hap2->alt_idxs[ins]);
    EXPECT_EQ(hap2->alt_idxs[ins], hap2->alt_idxs[del]);
    EXPECT_EQ(REC_CPX, hap2->rec_idxs[ins]);
    EXPECT_EQ(hap2->rec_idxs[ins], hap2->rec_idxs[del]);
    EXPECT_EQ(2, hap2->ploidies[ins]);
    EXPECT_EQ(hap2->ploidies[ins], hap2->ploidies[del]);
}

// Callers with no source record (e.g. CIGAR-derived variants) get the unknown sentinels.
TEST(ProvenanceVectorDefaults, UnknownSentinels) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars(new ctgVariants("ctg1"));
    vars->add_var(var_fields{.pos = 10, .rlen = 1, .type = TYPE_SUB, .loc = BED_INSIDE, .ref = "A",
            .alt = "G", .orig_gt = GT_ALT1_ALT1, .gt_qual = 60, .var_qual = 60, .phase_set = 0});

    ASSERT_EQ(1, vars->n);
    EXPECT_EQ(-1, vars->rec_idxs[0]);
    EXPECT_EQ(-1, vars->alt_idxs[0]);
    EXPECT_EQ(0, vars->ploidies[0]);
}

/* get_vartype ************************************************************************************/

TEST(GetVartype, SubSnp) {
    GlobalsGuard guard;
    g.sv_threshold = 8;

    // a substitution is a SNP regardless of allele length, so the threshold is never consulted
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_SUB, std::string(20, 'A'),
            std::string(20, 'C'));
    EXPECT_EQ(VARTYPE_SNP, vars->get_vartype(0));
}

TEST(GetVartype, SmallInsIndel) {
    GlobalsGuard guard;
    g.sv_threshold = 8;
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_INS, "",
            std::string(g.sv_threshold-1, 'A'));
    EXPECT_EQ(VARTYPE_INDEL, vars->get_vartype(0));
}

TEST(GetVartype, BoundaryInsSv) {
    GlobalsGuard guard;
    g.sv_threshold = 8;

    // the comparison is a strict <, so an alt of exactly sv_threshold bases is an SV
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_INS, "",
            std::string(g.sv_threshold, 'A'));
    EXPECT_EQ(VARTYPE_SV, vars->get_vartype(0));
}

TEST(GetVartype, SmallDelIndel) {
    GlobalsGuard guard;
    g.sv_threshold = 8;
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_DEL,
            std::string(g.sv_threshold-1, 'A'), "");
    EXPECT_EQ(VARTYPE_INDEL, vars->get_vartype(0));
}

TEST(GetVartype, BoundaryDelSv) {
    GlobalsGuard guard;
    g.sv_threshold = 8;
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_DEL,
            std::string(g.sv_threshold, 'A'), "");
    EXPECT_EQ(VARTYPE_SV, vars->get_vartype(0));
}

TEST(GetVartype, InsUsesAlt) {
    GlobalsGuard guard;
    g.sv_threshold = 8;

    // an INS is sized by its alt, so a long ref does not promote it to an SV
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_INS, std::string(20, 'A'), "CC");
    EXPECT_EQ(VARTYPE_INDEL, vars->get_vartype(0));
}

TEST(GetVartype, DelUsesRef) {
    GlobalsGuard guard;
    g.sv_threshold = 8;

    // a DEL is sized by its ref, so a short alt does not demote it to an INDEL
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_DEL, std::string(20, 'A'), "CC");
    EXPECT_EQ(VARTYPE_SV, vars->get_vartype(0));
}

TEST(GetVartype, CpxFallsToSv) {
    GlobalsGuard guard;
    g.sv_threshold = 8;

    // TYPE_CPX matches neither the SUB nor the INS/DEL branch, so it falls through to SV
    std::shared_ptr<ctgVariants> vars = make_typed_var(TYPE_CPX, "AC", "GT");
    EXPECT_EQ(VARTYPE_SV, vars->get_vartype(0));
}

/* set_allele_errtype *****************************************************************************/

// On a query record the truth allele count comes from matched_gts, recovered by alignment, and the
// query allele count from orig_gts, the record's own call.

TEST(SetAlleleErrtype, QueryZeroToOne) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_REF);
    EXPECT_EQ(AC_ERR_0_TO_1, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_0_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, QueryZeroToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EQ(AC_ERR_0_TO_2, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_0_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, QueryOneToZero) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_REF_ALT1);
    EXPECT_EQ(AC_ERR_1_TO_0, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_1_TO_0, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, QueryOneToOne) {
    GlobalsGuard guard;

    // a heterozygous call on the opposite haplotype is still one allele called for one expected
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_REF, GT_REF_ALT1);
    EXPECT_EQ(AC_ERR_1_TO_1, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_1_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, QueryOneToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    EXPECT_EQ(AC_ERR_1_TO_2, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_1_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, QueryTwoToZero) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_0, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_2_TO_0, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, QueryTwoToOne) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_1, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_2_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, QueryTwoToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_2, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_2_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, RefrefRefrefUnknown) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_REF_REF);

    // seed a different value so that a return-without-store would leave it behind
    vars->ac_errtype[0] = AC_ERR_2_TO_2;

    EXPECT_EQ(AC_UNKNOWN, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_UNKNOWN, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, HaploidUnknown) {
    GlobalsGuard guard;

    // a haploid genotype is not a diploid allele count, so it falls through to AC_UNKNOWN
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1, GT_ALT1);
    vars->ac_errtype[0] = AC_ERR_1_TO_1;

    EXPECT_EQ(AC_UNKNOWN, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_UNKNOWN, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, UnparseableGenotypeIsUnknownNotZeroAlleles) {
    GlobalsGuard guard;

    // a genotype carrying no diploid allele count is unknown, never silently zero alleles
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_MISSING, GT_ALT1_REF);
    vars->ac_errtype[0] = AC_ERR_1_TO_1;

    EXPECT_EQ(AC_UNKNOWN, vars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_UNKNOWN, vars->ac_errtype[0]);
}

// On a truth record the two sides swap: the truth allele count is the record's own orig_gts and the
// query allele count is matched_gts, recovered by alignment. The value keeps its absolute
// truth-then-query direction, so a truth record can reach *_TO_0 but never 0_TO_*.

TEST(SetAlleleErrtype, TruthOneToZero) {
    GlobalsGuard guard;

    // one truth allele, matched by no query allele: a pure false negative
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_REF);
    EXPECT_EQ(AC_ERR_1_TO_0, vars->set_allele_errtype(0, false));
    EXPECT_EQ(AC_ERR_1_TO_0, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TruthOneToOne) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_REF);
    EXPECT_EQ(AC_ERR_1_TO_1, vars->set_allele_errtype(0, false));
    EXPECT_EQ(AC_ERR_1_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TruthOneToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_1_TO_2, vars->set_allele_errtype(0, false));
    EXPECT_EQ(AC_ERR_1_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TruthTwoToZero) {
    GlobalsGuard guard;

    // two truth alleles, matched by no query allele: a pure false negative
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EQ(AC_ERR_2_TO_0, vars->set_allele_errtype(0, false));
    EXPECT_EQ(AC_ERR_2_TO_0, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TruthTwoToOne) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_1, vars->set_allele_errtype(0, false));
    EXPECT_EQ(AC_ERR_2_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TruthTwoToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_2, vars->set_allele_errtype(0, false));
    EXPECT_EQ(AC_ERR_2_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, FalseHomozygousAgreesAcrossCallsets) {
    GlobalsGuard guard;

    // one truth allele called as two: both records report it, so GE reads '+' on both samples
    std::shared_ptr<ctgVariants> qvars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    std::shared_ptr<ctgVariants> tvars = make_gt_var(GT_ALT1_REF, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_1_TO_2, qvars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_1_TO_2, tvars->set_allele_errtype(0, false));
    EXPECT_EQ("+", ac_strs[tvars->ac_errtype[0]]);
}

TEST(SetAlleleErrtype, FalseHeterozygousAgreesAcrossCallsets) {
    GlobalsGuard guard;

    // two truth alleles called as one: both records report it, so GE reads '-' on both samples
    std::shared_ptr<ctgVariants> qvars = make_gt_var(GT_ALT1_REF, GT_ALT1_ALT1);
    std::shared_ptr<ctgVariants> tvars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    EXPECT_EQ(AC_ERR_2_TO_1, qvars->set_allele_errtype(0, true));
    EXPECT_EQ(AC_ERR_2_TO_1, tvars->set_allele_errtype(0, false));
    EXPECT_EQ("-", ac_strs[tvars->ac_errtype[0]]);
}

TEST(SetAlleleErrtype, TruthHaploidUnknown) {
    GlobalsGuard guard;

    // a haploid genotype is not a diploid allele count on either side
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1, GT_ALT1);
    vars->ac_errtype[0] = AC_ERR_1_TO_1;

    EXPECT_EQ(AC_UNKNOWN, vars->set_allele_errtype(0, false));
    EXPECT_EQ(AC_UNKNOWN, vars->ac_errtype[0]);
}

/* matched_gt_is_swapped **************************************************************************/

TEST(CalcgtIsSwapped, EqualFalse) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_ALT1);
    EXPECT_FALSE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, OrigHomFalse) {
    GlobalsGuard guard;

    // orig 1|1 reports data from both haplotypes, so haplotype order does not matter
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    EXPECT_FALSE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, OrigRefFalse) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_REF_ALT1);
    EXPECT_FALSE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, CalcRefFalse) {
    GlobalsGuard guard;

    // matched 0|0 reports no data, so haplotype order does not matter
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_REF);
    EXPECT_FALSE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, HetOpposite0110True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_REF);
    EXPECT_TRUE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, HetOpposite1001True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_REF, GT_REF_ALT1);
    EXPECT_TRUE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, Orig01CalcHomCreditHap1True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.9;
    vars->credit[HAP2][0] = 0.1;

    // orig 0|1 expects the allele on HAP2, so the better HAP1 credit is reported by swapping
    EXPECT_TRUE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, Orig01CalcHomCreditHap2False) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.1;
    vars->credit[HAP2][0] = 0.9;
    EXPECT_FALSE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, Orig10CalcHomCreditHap2True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_REF, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.1;
    vars->credit[HAP2][0] = 0.9;
    EXPECT_TRUE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, CreditTieFalse) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.5;
    vars->credit[HAP2][0] = 0.5;

    // the comparison is a strict >, so equal credits keep the original haplotype order
    EXPECT_FALSE(vars->matched_gt_is_swapped(0));
}

TEST(CalcgtIsSwapped, UnexpectedErrors) {
    GlobalsGuard guard;

    // a haploid orig_gt reaches no branch, so the final else reports an unexpected pair
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1, GT_REF_ALT1);
    EXPECT_EXIT(vars->matched_gt_is_swapped(0), testing::ExitedWithCode(1),
            "Unexpected orig/matched genotypes for variant");
}

/* var_on_hap *************************************************************************************/

TEST(VarOnHap, HaploidAltBoth) {
    GlobalsGuard guard;

    // a haploid alternate allele is reported on both haplotypes
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1, GT_REF_REF);
    EXPECT_TRUE(vars->var_on_hap(0, HAP1));
    EXPECT_TRUE(vars->var_on_hap(0, HAP2));
}

TEST(VarOnHap, Gt10) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_REF, GT_REF_REF);
    EXPECT_TRUE(vars->var_on_hap(0, HAP1));
    EXPECT_FALSE(vars->var_on_hap(0, HAP2));
}

TEST(VarOnHap, Gt01) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_REF);
    EXPECT_FALSE(vars->var_on_hap(0, HAP1));
    EXPECT_TRUE(vars->var_on_hap(0, HAP2));
}

TEST(VarOnHap, HomBoth) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_TRUE(vars->var_on_hap(0, HAP1));
    EXPECT_TRUE(vars->var_on_hap(0, HAP2));
}

TEST(VarOnHap, RefrefNone) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_REF_REF);
    EXPECT_FALSE(vars->var_on_hap(0, HAP1));
    EXPECT_FALSE(vars->var_on_hap(0, HAP2));
}

TEST(VarOnHap, HaploidRefNone) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF, GT_REF_REF);
    EXPECT_FALSE(vars->var_on_hap(0, HAP1));
    EXPECT_FALSE(vars->var_on_hap(0, HAP2));
}

TEST(VarOnHap, CalcFlagSelects) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_REF);

    EXPECT_FALSE(vars->var_on_hap(0, HAP1, false));
    EXPECT_TRUE(vars->var_on_hap(0, HAP2, false));
    EXPECT_TRUE(vars->var_on_hap(0, HAP1, true));
    EXPECT_FALSE(vars->var_on_hap(0, HAP2, true));
}

/* set_var_matched_gt_on_hap **********************************************************************/

TEST(SetVarCalcgtOnHap, RefrefSetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    vars->set_var_matched_gt_on_hap(0, HAP1, true);
    EXPECT_EQ(GT_ALT1_REF, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, RefrefSetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    vars->set_var_matched_gt_on_hap(0, HAP2, true);
    EXPECT_EQ(GT_REF_ALT1, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, RefrefUnsetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP1, false), testing::ExitedWithCode(1),
            "Variant matched_gt already unset");
}

TEST(SetVarCalcgtOnHap, RefrefUnsetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP2, false), testing::ExitedWithCode(1),
            "Variant matched_gt already unset");
}

TEST(SetVarCalcgtOnHap, RefaltSetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    vars->set_var_matched_gt_on_hap(0, HAP1, true);
    EXPECT_EQ(GT_ALT1_ALT1, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, RefaltSetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP2, true), testing::ExitedWithCode(1),
            "Variant matched_gt already set");
}

TEST(SetVarCalcgtOnHap, RefaltUnsetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP1, false), testing::ExitedWithCode(1),
            "Variant matched_gt already unset");
}

TEST(SetVarCalcgtOnHap, RefaltUnsetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    vars->set_var_matched_gt_on_hap(0, HAP2, false);
    EXPECT_EQ(GT_REF_REF, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltrefSetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Variant matched_gt already set");
}

TEST(SetVarCalcgtOnHap, AltrefSetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    vars->set_var_matched_gt_on_hap(0, HAP2, true);
    EXPECT_EQ(GT_ALT1_ALT1, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltrefUnsetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    vars->set_var_matched_gt_on_hap(0, HAP1, false);
    EXPECT_EQ(GT_REF_REF, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltrefUnsetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP2, false), testing::ExitedWithCode(1),
            "Variant matched_gt already unset");
}

TEST(SetVarCalcgtOnHap, AltaltSetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Variant matched_gt already set");
}

TEST(SetVarCalcgtOnHap, AltaltSetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP2, true), testing::ExitedWithCode(1),
            "Variant matched_gt already set");
}

TEST(SetVarCalcgtOnHap, AltaltUnsetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);

    // clearing HAP1 leaves the alternate allele on HAP2 only
    vars->set_var_matched_gt_on_hap(0, HAP1, false);
    EXPECT_EQ(GT_REF_ALT1, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltaltUnsetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    vars->set_var_matched_gt_on_hap(0, HAP2, false);
    EXPECT_EQ(GT_ALT1_REF, vars->matched_gts[0]);
}

TEST(SetVarCalcgtOnHap, MissingErrors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_MISSING);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Unexpected matched_gts value");
}

TEST(SetVarCalcgtOnHap, HaploidErrors) {
    GlobalsGuard guard;

    // a haploid matched_gt matches none of the four diploid states, so the default branch errors
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Unexpected matched_gts value");
}

TEST(SetVarCalcgtOnHap, ErrorsSuppressedWithIgnore) {
    GlobalsGuard guard;

    // every invalid (state, hap, set) triple from the four diploid states
    struct transition { gt_t matched_gt; hap_t hap; bool set; };
    const std::vector<transition> invalid = {
        {GT_REF_REF,   HAP1, false}, {GT_REF_REF,   HAP2, false},
        {GT_REF_ALT1,  HAP2, true},  {GT_REF_ALT1,  HAP1, false},
        {GT_ALT1_REF,  HAP1, true},  {GT_ALT1_REF,  HAP2, false},
        {GT_ALT1_ALT1, HAP1, true},  {GT_ALT1_ALT1, HAP2, true},
    };

    for (const transition & t : invalid) {
        std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, t.matched_gt);
        vars->set_var_matched_gt_on_hap(0, t.hap, t.set, true);
        EXPECT_EQ(t.matched_gt, vars->matched_gts[0])
                << "matched_gt " << int(t.matched_gt) << " hap " << idx(t.hap)
                << " set " << t.set;
    }
}

TEST(SetVarCalcgtOnHap, IgnoreDoesNotSuppressMissing) {
    GlobalsGuard guard;

    // ignore_errors guards only the four diploid states; the default branch always errors
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_MISSING);
    EXPECT_EXIT(vars->set_var_matched_gt_on_hap(0, HAP1, true, true), testing::ExitedWithCode(1),
            "Unexpected matched_gts value");
}

/* variantData constructor ************************************************************************/

TEST(VariantDataCtor, DefaultCallsetQuery) {
    GlobalsGuard guard;
    variantData vcf;
    EXPECT_EQ(QUERY, vcf.callset);
}

TEST(VariantDataCtor, TwoHapMaps) {
    GlobalsGuard guard;
    variantData vcf;
    ASSERT_EQ(size_t(HAPS), vcf.variants.size());
    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
        EXPECT_TRUE(vcf.variants[hap].empty()) << "hap " << idx(hap);
    }
}

TEST(VariantDataCtor, EmptyMembers) {
    GlobalsGuard guard;
    variantData vcf;
    EXPECT_EQ(nullptr, vcf.ref);
    EXPECT_EQ("", vcf.filename);
    EXPECT_EQ("", vcf.sample);
    EXPECT_TRUE(vcf.contigs.empty());
    EXPECT_TRUE(vcf.lengths.empty());
    EXPECT_TRUE(vcf.observed_ploidies.empty());
}

/* parse-time filtering, counters, and summary warnings *******************************************/

/**
 * @class ParseVariants
 * @brief Restores parse-relevant global settings to their defaults before each test.
 */
class ParseVariants : public testing::Test {
protected:
    /** @brief Sets the globals parse_variants reads, leaving the histogram printed. */
    void SetUp() override {
        g.verbosity = 1; // print the genotype histogram, suppress per-variant warnings
        g.bed_exists = false;
        g.min_qual = 0;
        g.max_size = 1000;
        g.filters.clear();
        g.filter_ids.clear();
    }

    GlobalsGuard guard; ///< Saves global state on construction and restores it on destruction
    TempDir dir;        ///< Owns each test's fixture VCF, log, and output VCF
};

/* header and validation **************************************************************************/

// Contig lengths are copied into the output VCF header, so a contig line without one is fatal.
// htslib supplies IDX itself, so only a missing length can trip this check.
TEST_F(ParseVariants, ContigLineWithoutLengthErrors) {
    vcf_opts opts = make_vcf_opts();
    opts.contigs = {"##contig=<ID=chr1>"};
    EXPECT_EXIT(parse_unredirected(dir, {record(100, "A", "G", "1|0")}, opts),
            testing::ExitedWithCode(1), "header contig line didn't have 'IDX' and 'length'");
}

TEST_F(ParseVariants, MultipleSamplesErrors) {
    vcf_opts opts = make_vcf_opts();
    opts.sample = "QUERY1\tQUERY2";
    EXPECT_EXIT(parse_unredirected(dir, {fmt_record(100, "A", "G", "GT:PS", "1|0:1\t0|1:1")},
            opts),
            testing::ExitedWithCode(1), "Expected 1 sample but found 2");
}

// A selected filter absent from the header keeps its -1 sentinel, which no record's filter can
// match, so warning about it is the only notice that every variant was then dropped.
TEST_F(ParseVariants, SelectedFilterAbsentWarns) {
    g.filters = {"LOWQ"};
    g.filter_ids = {-1};
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0")});
    EXPECT_TRUE(logged(r, "Filter 'LOWQ' not found in QUERY VCF"));
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, "1 variants failed FILTER in QUERY VCF, skipped"));
}

// Returning to a contig that was already left behind means the VCF is not sorted by contig. Left
// unrejected, the contig list would gain a second entry for it and write_summary_vcf() would then
// walk that contig twice, emitting each of its variants twice.
TEST_F(ParseVariants, UnsortedContigErrors) {
    vcf_opts opts = make_vcf_opts(QUERY, {"chr1", "chr2"});
    EXPECT_EXIT(parse_unredirected(dir, {record(100, "A", "G", "1|0"),
            record(100, "A", "G", "1|0", "chr2"), record(200, "A", "G", "1|0")}, opts),
            testing::ExitedWithCode(1), "contig 'chr1' already parsed");
}

// Interleaving is what the guard rejects, not the contig order itself: a file whose contigs appear
// in an order the header does not use is still sorted, so long as each contig's records are grouped.
TEST_F(ParseVariants, ContigsOutOfHeaderOrderParse) {
    vcf_opts opts = make_vcf_opts(QUERY, {"chr1", "chr2"});
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0", "chr2"),
            record(100, "A", "G", "1|0"), record(200, "A", "G", "1|0")}, opts);
    EXPECT_FALSE(logged(r, "already parsed"));
    EXPECT_EQ(std::vector<std::string>({"chr2", "chr1"}), r.vars->contigs);
    EXPECT_EQ(size_t(1), count_pos(r, 100));
    EXPECT_EQ(size_t(1), count_pos(r, 200));
}

// Documents rather than enforces: the seqnames failure at variant.cpp:829 is unreachable. An empty
// contig dictionary still yields a non-NULL array, and ERROR() exits, so its `goto error1` is dead.
TEST_F(ParseVariants, HeaderWithoutContigsParsesNoContigs) {
    vcf_opts opts = make_vcf_opts();
    opts.contigs = {};
    ParseResult r = parse_records(dir, {}, opts);
    EXPECT_TRUE(r.vars->contigs.empty());
    EXPECT_TRUE(r.vars->variants[HAP1].empty());
    EXPECT_FALSE(logged(r, "Failed to read QUERY VCF"));
}

/* filter and quality *****************************************************************************/

// A record whose only filter is one the user did not select is dropped, and one that PASSes is not.
TEST_F(ParseVariants, FilterFailSkipped) {
    g.filters = {"PASS"};
    g.filter_ids = {-1}; // parse_variants fills this in from the header's FILTER IDX
    vcf_opts opts = make_vcf_opts();
    opts.filters = {"##FILTER=<ID=PASS,Description=\"All filters passed\">",
                    "##FILTER=<ID=LOWQ,Description=\"Low quality\">"};
    ParseResult r = parse_records(dir, {qual_filter_record(100, "50", "LOWQ"),
                                        qual_filter_record(200, "50", "PASS")}, opts);
    EXPECT_EQ(1, total_kept(r));
    EXPECT_FALSE(kept_pos(r, 100));
    EXPECT_TRUE(kept_pos(r, 200));
    EXPECT_TRUE(logged(r, "1 variants failed FILTER in QUERY VCF, skipped"));
}

// With no filters selected, the FILTER column is not consulted at all.
TEST_F(ParseVariants, NoFiltersSelectedPasses) {
    vcf_opts opts = make_vcf_opts();
    opts.filters = {"##FILTER=<ID=PASS,Description=\"All filters passed\">",
                    "##FILTER=<ID=LOWQ,Description=\"Low quality\">"};
    ParseResult r = parse_records(dir, {qual_filter_record(100, "50", "LOWQ")}, opts);
    EXPECT_EQ(1, total_kept(r));
    EXPECT_FALSE(logged(r, "failed FILTER"));
}

// A record with no filters at all passes even when a filter was selected, since there is nothing
// to contradict the selection.
TEST_F(ParseVariants, UnfilteredRecordPassesSelectedFilter) {
    g.filters = {"PASS"};
    g.filter_ids = {-1};
    ParseResult r = parse_records(dir, {qual_filter_record(100, "50", ".")});
    EXPECT_EQ(1, total_kept(r));
    EXPECT_FALSE(logged(r, "failed FILTER"));
}

// The comparison is inclusive, so a variant exactly at min_qual is kept and one below is not.
TEST_F(ParseVariants, BelowMinQualSkipped) {
    g.min_qual = 60;
    ParseResult r = parse_records(dir, {qual_filter_record(100, "59", "PASS"),
                                        qual_filter_record(200, "60", "PASS")});
    EXPECT_EQ(1, total_kept(r));
    EXPECT_FALSE(kept_pos(r, 100));
    EXPECT_TRUE(kept_pos(r, 200));
    EXPECT_TRUE(logged(r, "1 variants of low quality (<60) in QUERY VCF, skipped"));
}

// An unreported QUAL parses as NaN, which no comparison would accept, so it is read as zero.
TEST_F(ParseVariants, QualNanReadAsZero) {
    ParseResult r = parse_records(dir, {qual_filter_record(100, ".", "PASS")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_FLOAT_EQ(0, hap_vars(r, HAP1)->var_quals[0]);
}

TEST_F(ParseVariants, QualNanSkippedWhenMinQualPositive) {
    g.min_qual = 1;
    ParseResult r = parse_records(dir, {qual_filter_record(100, ".", "PASS")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, "1 variants of low quality (<1) in QUERY VCF, skipped"));
}

// An integer GQ becomes the variant's genotype quality, which is stored unclamped.
TEST_F(ParseVariants, IntegerGqStored) {
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT:GQ", "1|0:44")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_FLOAT_EQ(44, hap_vars(r, HAP1)->gt_quals[0]);
}

// A float GQ makes the integer read fail with -2, and the float retry truncates toward zero.
TEST_F(ParseVariants, FloatGqTruncatedToInt) {
    vcf_opts opts = make_vcf_opts();
    opts.formats = {"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                    "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality\">"};
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT:GQ", "1|0:33.7")}, opts);
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_FLOAT_EQ(33, hap_vars(r, HAP1)->gt_quals[0]);
}

// A GQ declared in the header but absent from a record is read as zero, not as an error.
TEST_F(ParseVariants, MissingGqReadAsZero) {
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_FLOAT_EQ(0, hap_vars(r, HAP1)->gt_quals[0]);
}

/* genotype ***************************************************************************************/

// Without a GT declaration there is no genotype to read, so every record is one haploid alternate.
TEST_F(ParseVariants, NoGtInHeaderWarnsAndAssumesMonoploid) {
    vcf_opts opts = make_vcf_opts();
    opts.formats = {"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">"};
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GQ:PS", "44:1")}, opts);
    EXPECT_TRUE(logged(r, "'GT' tag not defined in QUERY VCF header, assuming monoploid"));
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(0, hap_vars(r, HAP2)->n);
    EXPECT_EQ(1, hap_vars(r, HAP1)->ploidies[0]);
    EXPECT_EQ(GT_ALT1_REF, hap_vars(r, HAP1)->orig_gts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_ALT1, 1)));
    EXPECT_EQ(std::vector< std::set<int> >({{1}}), r.vars->observed_ploidies);
}

// A haploid alternate call is reported on haplotype 1 alone, with ploidy 1 recorded.
TEST_F(ParseVariants, HaploidAltKeptOnHap1) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(0, hap_vars(r, HAP2)->n);
    EXPECT_EQ(GT_ALT1_REF, hap_vars(r, HAP1)->orig_gts[0]);
    EXPECT_EQ(1, hap_vars(r, HAP1)->ploidies[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_ALT1, 1)));
}

TEST_F(ParseVariants, HaploidRefSkipped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, gt_hist_line(GT_REF, 1)));
}

TEST_F(ParseVariants, Gt01KeptOnHap2) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "0|1")});
    EXPECT_EQ(0, hap_vars(r, HAP1)->n);
    ASSERT_EQ(1, hap_vars(r, HAP2)->n);
    EXPECT_EQ(GT_REF_ALT1, hap_vars(r, HAP2)->orig_gts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_REF_ALT1, 1)));
}

TEST_F(ParseVariants, Gt10KeptOnHap1) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(0, hap_vars(r, HAP2)->n);
    EXPECT_EQ(GT_ALT1_REF, hap_vars(r, HAP1)->orig_gts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_ALT1_REF, 1)));
}

// A homozygous record yields the same allele on both haplotypes.
TEST_F(ParseVariants, Gt11SplitAcrossBothHaps) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|1")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    ASSERT_EQ(1, hap_vars(r, HAP2)->n);
    EXPECT_EQ("G", hap_vars(r, HAP1)->alts[0]);
    EXPECT_EQ("G", hap_vars(r, HAP2)->alts[0]);
    EXPECT_EQ(99, hap_vars(r, HAP1)->poss[0]);
    EXPECT_EQ(99, hap_vars(r, HAP2)->poss[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_ALT1_ALT1, 1)));
    EXPECT_TRUE(logged(r, "1 homozygous and multi-allelic variants in QUERY VCF, split"));
}

// A homozygous record contributes one variant per haplotype, and each keeps the 1|1 label unless
// the other haplotype's copy was dropped. With one record and nothing to overlap, both keep it.
TEST_F(ParseVariants, HomozygousKeepsBothAllelesWhenNothingSkipped) {
    const std::vector< std::pair<std::string, std::string> > alleles = {
        {"A", "G"},    // SUB, rlen 1
        {"AGG", "A"},  // DEL, rlen 2
        {"A", "AGG"},  // INS, rlen 0
    };
    for (const auto & [ref, alt] : alleles) {
        ParseResult r = parse_records(dir, {record(100, ref, alt, "1|1")});
        ASSERT_EQ(1, hap_vars(r, HAP1)->n) << ref << " -> " << alt;
        ASSERT_EQ(1, hap_vars(r, HAP2)->n) << ref << " -> " << alt;
        EXPECT_EQ(GT_ALT1_ALT1, hap_vars(r, HAP1)->orig_gts[0]) << ref << " -> " << alt;
        EXPECT_EQ(GT_ALT1_ALT1, hap_vars(r, HAP2)->orig_gts[0]) << ref << " -> " << alt;
        EXPECT_FALSE(logged(r, "overlapping variants")) << ref << " -> " << alt;
    }
}

// A compound heterozygote puts a different ALT on each haplotype, each relabelled as a simple het.
TEST_F(ParseVariants, Gt12SplitByAllele) {
    ParseResult r = parse_records(dir, {record(100, "A", "G,T", "1|2")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    ASSERT_EQ(1, hap_vars(r, HAP2)->n);
    EXPECT_EQ("G", hap_vars(r, HAP1)->alts[0]);
    EXPECT_EQ("T", hap_vars(r, HAP2)->alts[0]);
    EXPECT_EQ(GT_ALT1_REF, hap_vars(r, HAP1)->orig_gts[0]);
    EXPECT_EQ(GT_REF_ALT1, hap_vars(r, HAP2)->orig_gts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_ALT1_ALT2, 1)));
}

TEST_F(ParseVariants, Gt21SplitByAllele) {
    ParseResult r = parse_records(dir, {record(100, "A", "G,T", "2|1")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    ASSERT_EQ(1, hap_vars(r, HAP2)->n);
    EXPECT_EQ("T", hap_vars(r, HAP1)->alts[0]);
    EXPECT_EQ("G", hap_vars(r, HAP2)->alts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_ALT2_ALT1, 1)));
}

// 0|2 reaches no named genotype, so it is tallied as other while its second ALT is still kept.
TEST_F(ParseVariants, Gt02TalliedAsOther) {
    ParseResult r = parse_records(dir, {record(100, "A", "G,T", "0|2")});
    EXPECT_EQ(0, hap_vars(r, HAP1)->n);
    ASSERT_EQ(1, hap_vars(r, HAP2)->n);
    EXPECT_EQ("T", hap_vars(r, HAP2)->alts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_OTHER, 1)));
}

TEST_F(ParseVariants, PolyploidErrors) {
    EXPECT_EXIT(parse_unredirected(dir, {record(100, "A", "G", "1|1|1")}, make_vcf_opts()),
            testing::ExitedWithCode(1), "found variant with ploidy 3");
}

// Every record's ploidy is recorded, and a contig carrying more than one is not an error.
TEST_F(ParseVariants, MixedPloidyRecordedAndKept) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|1"),
                                        record(200, "A", "G", "1")});
    EXPECT_EQ(std::vector< std::set<int> >({{1, 2}}), r.vars->observed_ploidies);
    EXPECT_FALSE(logged(r, "ploidy"));
    EXPECT_EQ(2, hap_vars(r, HAP1)->n);
    EXPECT_EQ(1, hap_vars(r, HAP2)->n);
}

// chrX was previously exempted from the mismatch warning by name; nothing warns now, so the
// exemption is gone and the contig is treated no differently from any other.
TEST_F(ParseVariants, MixedPloidyOnChrXSilent) {
    vcf_opts opts = make_vcf_opts(QUERY, {"chrX"});
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|1", "chrX"),
                                        record(200, "A", "G", "1", "chrX")}, opts);
    EXPECT_FALSE(logged(r, "ploidy"));
    EXPECT_EQ(2, r.vars->variants[HAP1]["chrX"]->n);
}

// chrY never had that exemption, so it warned on every record after the first. It no longer does.
TEST_F(ParseVariants, MixedPloidyOnChrYSilent) {
    vcf_opts opts = make_vcf_opts(QUERY, {"chrY"});
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|1", "chrY"),
                                        record(200, "A", "G", "1", "chrY")}, opts);
    EXPECT_FALSE(logged(r, "ploidy"));
    EXPECT_EQ(2, r.vars->variants[HAP1]["chrY"]->n);
}

/* missing (.) alleles ****************************************************************************/

// A no-call has no known allele on either haplotype, so the whole record is dropped.
TEST_F(ParseVariants, NoCallDroppedAndCounted) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", ".|.")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(kept_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants with no known alleles (.|.) in QUERY VCF, skipped"));
    EXPECT_TRUE(logged(r, gt_hist_line(GT_MISSING, 1)));
}

// A no-call is one dropped record, not one dropped record per haplotype.
TEST_F(ParseVariants, NoCallCountedOncePerRecord) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", ".|.")});
    EXPECT_TRUE(logged(r, "1 variants with no known alleles"));
    EXPECT_FALSE(logged(r, "2 variants with no known alleles"));
}

// A half call keeps its known allele, so it must not be tallied as a no-call.
TEST_F(ParseVariants, HalfCallCountedDistinctlyFromNoCall) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|.")});
    EXPECT_TRUE(logged(r, gt_hist_line(GT_HALF, 1)));
    EXPECT_FALSE(logged(r, ".|.:")); // histogram line is only printed for nonzero counts
    EXPECT_FALSE(logged(r, "no known alleles"));
}

// The summary must not claim a half call was skipped, because its known allele was evaluated.
TEST_F(ParseVariants, HalfCallNotReportedAsSkipped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|.")});
    EXPECT_EQ(1, total_kept(r));
    EXPECT_TRUE(kept_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants with a half call (1|.) in QUERY VCF, known allele kept"));
    EXPECT_FALSE(logged(r, "skipped"));
}

// A missing allele on either haplotype leaves the known allele on the other.
TEST_F(ParseVariants, HalfCallKeptOnTheHaplotypeWithTheKnownAllele) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|."),
                                        record(200, "A", "G", ".|1")});
    EXPECT_EQ(1, kept_on_hap(r, HAP1));
    EXPECT_EQ(1, kept_on_hap(r, HAP2));
    EXPECT_TRUE(logged(r, gt_hist_line(GT_HALF, 2)));
    EXPECT_TRUE(logged(r, "2 variants with a half call"));
}

// A half call keeps its known allele on the haplotype that carries it; the phase bit htslib sets on
// a phased '.' is what clears the unphased-heterozygote guard.
TEST_F(ParseVariants, HalfCallHap1AlleleKept) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|.")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(0, hap_vars(r, HAP2)->n);
    EXPECT_EQ(GT_ALT1_REF, hap_vars(r, HAP1)->orig_gts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_HALF, 1)));
}

TEST_F(ParseVariants, HalfCallHap2AlleleKept) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", ".|1")});
    EXPECT_EQ(0, hap_vars(r, HAP1)->n);
    ASSERT_EQ(1, hap_vars(r, HAP2)->n);
    EXPECT_EQ(GT_REF_ALT1, hap_vars(r, HAP2)->orig_gts[0]);
    EXPECT_TRUE(logged(r, gt_hist_line(GT_HALF, 1)));
}

// An unphased half call is still a half call, but its known allele is dropped as unphased.
TEST_F(ParseVariants, UnphasedHalfCallSkipped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "./1")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, gt_hist_line(GT_HALF, 1)));
    EXPECT_TRUE(logged(r, "1 variants with a half call (1|.) in QUERY VCF, known allele kept"));
    EXPECT_TRUE(logged(r, "1 variants with unphased genotypes in QUERY VCF, skipped"));
}

// A haploid no-call is a no-call, not an alternate: the missing-allele test precedes the truthiness
// test on bcf_gt_allele()'s -1, which would otherwise read as ALT.
TEST_F(ParseVariants, HaploidNoCallSkippedNotTalliedAsAlt) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", ".")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, gt_hist_line(GT_MISSING, 1)));
    EXPECT_FALSE(logged(r, gt_hist_line(GT_ALT1, 1)));
    EXPECT_TRUE(logged(r, "1 variants with no known alleles (.|.) in QUERY VCF, skipped"));
}

// The per-variant warning must say the known allele is kept, since only the missing one is dropped.
TEST_F(ParseVariants, HalfCallPerVariantWarningKeepsKnownAllele) {
    g.verbosity = 2; // per-variant warnings are only printed when verbose
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|.")});
    EXPECT_TRUE(logged(r, "Variant with a half call (1|.) in QUERY VCF at chr1:99,"
            " keeping known allele"));
}

TEST_F(ParseVariants, NoCallPerVariantWarningSaysSkipping) {
    g.verbosity = 2;
    ParseResult r = parse_records(dir, {record(100, "A", "G", ".|.")});
    EXPECT_TRUE(logged(r, "Variant with no known alleles (.|.) in QUERY VCF at chr1:99, skipping"));
}

// A spanning deletion is tallied as a reference-type allele, whereas a missing allele is tallied
// only as a genotype, so the two paths cannot be confused in the type summary.
TEST_F(ParseVariants, SpanningDeletionTalliedAsRefTypeButNoCallIsNot) {
    ParseResult star = parse_records(dir, {record(100, "A", "*", "1|0")});
    EXPECT_TRUE(logged(star, type_hist_line(TYPE_REF, 1)));

    ParseResult missing = parse_records(dir, {record(100, "A", "G", ".|.")});
    EXPECT_FALSE(logged(missing, type_hist_line(TYPE_REF, 1)));
    EXPECT_TRUE(logged(missing, gt_hist_line(GT_MISSING, 1)));
}

/* phase set **************************************************************************************/

// Without a PS declaration every variant on a contig shares one implicit phase set, numbered zero.
TEST_F(ParseVariants, PhaseSetNotInHeaderWarns) {
    vcf_opts opts = make_vcf_opts();
    opts.formats = {"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"};
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT", "1|0")}, opts);
    EXPECT_TRUE(logged(r,
            "'PS' tag not defined in QUERY VCF header, assuming one phase set per contig"));
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(0, hap_vars(r, HAP1)->phase_sets[0]);
}

// A heterozygote without a PS tag has an unknown phase set, so it is counted and reported.
TEST_F(ParseVariants, PhaseSetMissingOnHetWarns) {
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT", "1|0")});
    EXPECT_TRUE(logged(r, "1 variants missing PS tags in QUERY VCF, kept"));
    ASSERT_EQ(1, hap_vars(r, HAP1)->n); // counted, not dropped
    EXPECT_EQ(0, hap_vars(r, HAP1)->phase_sets[0]);
}

TEST_F(ParseVariants, PhaseSetPresentStored) {
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT:PS", "1|0:7")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(7, hap_vars(r, HAP1)->phase_sets[0]);
}

// A homozygote is exempt from the missing-PS count, since its phase needs no resolving, but it
// still declared no phase set of its own and must not be attributed to the previous record's.
TEST_F(ParseVariants, PhaseSetNotInheritedByHomRecordWithoutPs) {
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT:PS", "1|0:7"),
                                        fmt_record(200, "A", "G", "GT", "1|1")});
    ASSERT_EQ(2, hap_vars(r, HAP1)->n);
    EXPECT_EQ(7, hap_vars(r, HAP1)->phase_sets[0]);
    EXPECT_EQ(0, hap_vars(r, HAP1)->phase_sets[1]);
    EXPECT_FALSE(logged(r, "missing PS tags")); // a homozygote needs no phase set to be resolved
}

// The same holds for a haploid record, also exempt from the count and also owed no inherited value.
TEST_F(ParseVariants, PhaseSetNotInheritedByHaploidRecordWithoutPs) {
    ParseResult r = parse_records(dir, {fmt_record(100, "A", "G", "GT:PS", "1:7"),
                                        fmt_record(200, "A", "G", "GT", "1")});
    ASSERT_EQ(2, hap_vars(r, HAP1)->n);
    EXPECT_EQ(7, hap_vars(r, HAP1)->phase_sets[0]);
    EXPECT_EQ(0, hap_vars(r, HAP1)->phase_sets[1]);
    EXPECT_FALSE(logged(r, "missing PS tags")); // a haploid variant needs no phase set either
}

/* allele filtering *******************************************************************************/

// Unphased heterozygous genotypes are still dropped at parse time, counter and warning intact.
TEST_F(ParseVariants, UnphasedHeterozygousGenotypeStillDropped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "0/1")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(kept_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants with unphased genotypes in QUERY VCF, skipped"));
}

// A spanning deletion allele carries no variation, so it is dropped and counted.
TEST_F(ParseVariants, SpanningDeletionDroppedAndCounted) {
    ParseResult r = parse_records(dir, {record(100, "A", "*", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(kept_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants spanned by deletion in QUERY VCF, skipped"));
}

// An unphased homozygote states the same allele on both haplotypes, so its phase is not in doubt.
TEST_F(ParseVariants, UnphasedHomAllowed) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1/1")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    ASSERT_EQ(1, hap_vars(r, HAP2)->n);
    EXPECT_EQ(GT_ALT1_ALT1, hap_vars(r, HAP1)->orig_gts[0]);
    EXPECT_FALSE(logged(r, "unphased genotypes"));
}

/* classification: prefix/suffix trim and CPX split ***********************************************/

// An ALT identical to its REF carries no variation, whether one base long or several.
TEST_F(ParseVariants, RefCallDroppedAndCounted) {
    ParseResult r = parse_records(dir, {record(100, "A", "A", "1|0"),
                                        record(200, "AT", "AT", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(kept_pos(r, 100));
    EXPECT_FALSE(kept_pos(r, 200));
    EXPECT_TRUE(logged(r, "2 reference variants in QUERY VCF, skipped"));
}

TEST_F(ParseVariants, SnpClassified) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(TYPE_SUB, hap_vars(r, HAP1)->types[0]);
    EXPECT_EQ(99, hap_vars(r, HAP1)->poss[0]);
    EXPECT_EQ(1, hap_vars(r, HAP1)->rlens[0]);
    EXPECT_EQ("A", hap_vars(r, HAP1)->refs[0]);
    EXPECT_EQ("G", hap_vars(r, HAP1)->alts[0]);
    EXPECT_TRUE(logged(r, type_hist_line(TYPE_SUB, 1)));
}

// An insertion drops the anchor base it shares with REF and advances past it.
TEST_F(ParseVariants, InsertionTrimsSharedPrefix) {
    ParseResult r = parse_records(dir, {record(100, "A", "AGG", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(TYPE_INS, hap_vars(r, HAP1)->types[0]);
    EXPECT_EQ(100, hap_vars(r, HAP1)->poss[0]);
    EXPECT_EQ(0, hap_vars(r, HAP1)->rlens[0]);
    EXPECT_EQ("", hap_vars(r, HAP1)->refs[0]);
    EXPECT_EQ("GG", hap_vars(r, HAP1)->alts[0]);
}

TEST_F(ParseVariants, DeletionTrimsSharedPrefix) {
    ParseResult r = parse_records(dir, {record(100, "AGG", "A", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(TYPE_DEL, hap_vars(r, HAP1)->types[0]);
    EXPECT_EQ(100, hap_vars(r, HAP1)->poss[0]);
    EXPECT_EQ(2, hap_vars(r, HAP1)->rlens[0]);
    EXPECT_EQ("GG", hap_vars(r, HAP1)->refs[0]);
    EXPECT_EQ("", hap_vars(r, HAP1)->alts[0]);
}

// Equal-length alleles agreeing past their first base are one substitution, not a complex variant.
TEST_F(ParseVariants, MnpWithSharedSuffixIsSub) {
    ParseResult r = parse_records(dir, {record(100, "AT", "GT", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(TYPE_SUB, hap_vars(r, HAP1)->types[0]);
    EXPECT_EQ(99, hap_vars(r, HAP1)->poss[0]);
    EXPECT_EQ(1, hap_vars(r, HAP1)->rlens[0]);
    EXPECT_EQ("A", hap_vars(r, HAP1)->refs[0]);
    EXPECT_EQ("G", hap_vars(r, HAP1)->alts[0]);
}

// Equal-length alleles that share no suffix are complex, so they become a co-located INS and DEL.
TEST_F(ParseVariants, ComplexEqualLengthSplitIntoInsAndDel) {
    ParseResult r = parse_records(dir, {record(100, "AT", "GC", "1|0")});
    std::shared_ptr<ctgVariants> h1 = hap_vars(r, HAP1);
    ASSERT_EQ(2, h1->n);
    EXPECT_EQ(TYPE_INS, h1->types[0]);
    EXPECT_EQ(99, h1->poss[0]);
    EXPECT_EQ(0, h1->rlens[0]);
    EXPECT_EQ("", h1->refs[0]);
    EXPECT_EQ("GC", h1->alts[0]);
    EXPECT_EQ(TYPE_DEL, h1->types[1]);
    EXPECT_EQ(99, h1->poss[1]);
    EXPECT_EQ(2, h1->rlens[1]);
    EXPECT_EQ("AT", h1->refs[1]);
    EXPECT_EQ("", h1->alts[1]);
    EXPECT_TRUE(logged(r, "1 complex (CPX) variants in QUERY VCF, split into INS + DEL"));

    // the type summary counts the original CPX allele, not the two variants it became
    EXPECT_TRUE(logged(r, type_hist_line(TYPE_CPX, 1)));
}

// A net insertion whose flanks do not match is complex, and keeps its untrimmed alleles.
TEST_F(ParseVariants, ComplexInsertionSplitIntoInsAndDel) {
    ParseResult r = parse_records(dir, {record(100, "AT", "GCC", "1|0")});
    std::shared_ptr<ctgVariants> h1 = hap_vars(r, HAP1);
    ASSERT_EQ(2, h1->n);
    EXPECT_EQ(TYPE_INS, h1->types[0]);
    EXPECT_EQ("GCC", h1->alts[0]);
    EXPECT_EQ(TYPE_DEL, h1->types[1]);
    EXPECT_EQ("AT", h1->refs[1]);
    EXPECT_EQ(2, h1->rlens[1]);
    EXPECT_EQ(99, h1->poss[0]);
    EXPECT_EQ(99, h1->poss[1]);
}

TEST_F(ParseVariants, ComplexDeletionSplitIntoInsAndDel) {
    ParseResult r = parse_records(dir, {record(100, "ATG", "CC", "1|0")});
    std::shared_ptr<ctgVariants> h1 = hap_vars(r, HAP1);
    ASSERT_EQ(2, h1->n);
    EXPECT_EQ(TYPE_INS, h1->types[0]);
    EXPECT_EQ("CC", h1->alts[0]);
    EXPECT_EQ(TYPE_DEL, h1->types[1]);
    EXPECT_EQ("ATG", h1->refs[1]);
    EXPECT_EQ(3, h1->rlens[1]);
}

// Soft-masked reference sequence reaches the parser in lowercase, and is stored uppercased.
TEST_F(ParseVariants, LowercaseAllelesUppercased) {
    ParseResult r = parse_records(dir, {record(100, "a", "g", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ("A", hap_vars(r, HAP1)->refs[0]);
    EXPECT_EQ("G", hap_vars(r, HAP1)->alts[0]);
}

// Case carries no biological meaning, so an ALT matching its REF apart from case is a reference
// call, not a substitution of a base for itself.
TEST_F(ParseVariants, CaseOnlyRefCallDroppedAndCounted) {
    ParseResult r = parse_records(dir, {record(100, "a", "A", "1|0"),
                                        record(200, "A", "a", "1|0"),
                                        record(300, "at", "AT", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(kept_pos(r, 100));
    EXPECT_FALSE(kept_pos(r, 200));
    EXPECT_FALSE(kept_pos(r, 300));
    EXPECT_TRUE(logged(r, "3 reference variants in QUERY VCF, skipped"));
}

// A lowercase anchor base is still an anchor base, so it must trim and leave a plain insertion.
TEST_F(ParseVariants, LowercaseAnchorTrimsToInsertion) {
    ParseResult r = parse_records(dir, {record(100, "A", "aGG", "1|0")});
    std::shared_ptr<ctgVariants> h1 = hap_vars(r, HAP1);
    ASSERT_EQ(1, h1->n);
    EXPECT_EQ(TYPE_INS, h1->types[0]);
    EXPECT_EQ(100, h1->poss[0]); // 0-based 99, past the trimmed anchor base
    EXPECT_EQ("", h1->refs[0]);
    EXPECT_EQ("GG", h1->alts[0]);
    EXPECT_FALSE(logged(r, "complex (CPX) variants"));
}

// The same holds for a deletion, whose anchor base is the whole of a lowercase ALT.
TEST_F(ParseVariants, LowercaseAnchorTrimsToDeletion) {
    ParseResult r = parse_records(dir, {record(100, "ATT", "a", "1|0")});
    std::shared_ptr<ctgVariants> h1 = hap_vars(r, HAP1);
    ASSERT_EQ(1, h1->n);
    EXPECT_EQ(TYPE_DEL, h1->types[0]);
    EXPECT_EQ(100, h1->poss[0]); // 0-based 99, past the trimmed anchor base
    EXPECT_EQ("TT", h1->refs[0]);
    EXPECT_EQ("", h1->alts[0]);
    EXPECT_FALSE(logged(r, "complex (CPX) variants"));
}

// A shared suffix differing only in case must still chop, leaving a SNP rather than a CPX.
TEST_F(ParseVariants, LowercaseSuffixTrimsToSubstitution) {
    ParseResult r = parse_records(dir, {record(100, "At", "GT", "1|0")});
    std::shared_ptr<ctgVariants> h1 = hap_vars(r, HAP1);
    ASSERT_EQ(1, h1->n);
    EXPECT_EQ(TYPE_SUB, h1->types[0]);
    EXPECT_EQ(99, h1->poss[0]);
    EXPECT_EQ("A", h1->refs[0]);
    EXPECT_EQ("G", h1->alts[0]);
    EXPECT_FALSE(logged(r, "complex (CPX) variants"));
}

// Documents rather than enforces: the unexpected-type error at variant.cpp:1130 looks unreachable,
// because each of the three allele-length branches either assigns a type or skips the allele. Every
// shape below therefore parses without exiting, and stores only the three types add_var() can
// receive -- never TYPE_CPX, which is always split into an insertion and a deletion first.
TEST_F(ParseVariants, EveryAlleleShapeReachesAKnownType) {
    const std::vector< std::pair<std::string, std::string> > alleles = {
        {"A", "G"}, {"AT", "GT"}, {"AT", "GC"}, {"A", "AGG"}, {"AGG", "A"},
        {"AT", "GCC"}, {"ATG", "CC"}, {"A", "<DEL>"},
    };
    for (const auto & [ref, alt] : alleles) {
        ParseResult r = parse_records(dir, {record(100, ref, alt, "1|0")});
        ASSERT_LE(1, hap_vars(r, HAP1)->n) << ref << " -> " << alt;
        for (int vi = 0; vi < hap_vars(r, HAP1)->n; vi++) {
            const edittype_t type = hap_vars(r, HAP1)->types[vi];
            EXPECT_TRUE(type == TYPE_SUB || type == TYPE_INS || type == TYPE_DEL)
                    << ref << " -> " << alt << " variant " << vi << " type " << int(type);
        }
    }
}

/* region, size, and overlap **********************************************************************/

TEST_F(ParseVariants, InsideRegionKept) {
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{50, 200}});
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0")});
    ASSERT_EQ(1, hap_vars(r, HAP1)->n);
    EXPECT_EQ(BED_INSIDE, hap_vars(r, HAP1)->locs[0]);
}

TEST_F(ParseVariants, OutsideRegionSkipped) {
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{200, 300}});
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, "1 variants outside selected regions in QUERY VCF, skipped"));
}

// A variant straddling a region boundary is only partly evaluable, so it is dropped as a border
// case.
TEST_F(ParseVariants, BorderRegionSkipped) {
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{100, 200}});
    ParseResult r = parse_records(dir, {record(100, "AGGG", "A", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, "1 variants on border of selected regions in QUERY VCF, skipped"));
}

// A contig absent from the BED was not selected at all, and is reported alongside outside variants.
TEST_F(ParseVariants, OffContigSkipped) {
    g.bed_exists = true;
    g.bed = make_bed("chr2", {{50, 200}});
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_TRUE(logged(r, "1 variants outside selected regions in QUERY VCF, skipped"));
}

// The size limit is applied to the trimmed alleles, so the anchor base does not count toward it.
TEST_F(ParseVariants, TooLargeVariantSkipped) {
    g.max_size = 5;
    ParseResult r = parse_records(dir, {record(100, "A", "A" + std::string(10, 'G'), "1|0"),
                                        record(200, "A", "A" + std::string(5, 'G'), "1|0")});
    EXPECT_EQ(1, total_kept(r));
    EXPECT_TRUE(logged(r, "1 large (size > 5) variants in QUERY VCF, skipped"));
}

// Overlapping variants are still dropped at parse time, with their counter and warning intact.
TEST_F(ParseVariants, OverlappingVariantStillDropped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0"),
                                        record(100, "A", "T", "1|0")});
    EXPECT_EQ(1, kept_on_hap(r, HAP1));
    EXPECT_EQ(0, kept_on_hap(r, HAP2));
    EXPECT_EQ("G", hap_vars(r, HAP1)->alts[0]); // the first ALT, not the overlapping "T"
    EXPECT_TRUE(logged(r, "1 overlapping variants in QUERY VCF, skipped"));
}

// Two insertions at one position would both anchor to the same base, so the second is dropped even
// though it consumes no reference and therefore does not overlap the first.
TEST_F(ParseVariants, TwoInsertionsAtSamePositionSkipped) {
    ParseResult r = parse_records(dir, {record(100, "A", "AGG", "1|0"),
                                        record(100, "A", "ATT", "1|0")});
    EXPECT_EQ(1, kept_on_hap(r, HAP1));
    EXPECT_EQ("GG", hap_vars(r, HAP1)->alts[0]);
    EXPECT_TRUE(logged(r, "1 overlapping variants in QUERY VCF, skipped"));
}

// When one half of a homozygous record is dropped as overlapping, the surviving half stops claiming
// to be homozygous, so its allele is not reported on a haplotype that has no variant. Here the
// haplotype-1 copy overlaps a preceding deletion, leaving the haplotype-2 copy as 0|1.
TEST_F(ParseVariants, HomozygousDowngradedWhenOneHapOverlaps) {
    ParseResult r = parse_records(dir, {record(100, "AGGT", "A", "1|0"),
                                        record(102, "C", "T", "1|1")});
    ASSERT_EQ(1, kept_on_hap(r, HAP1)); // the deletion only; its SNP half overlapped
    ASSERT_EQ(1, kept_on_hap(r, HAP2));
    EXPECT_EQ(101, hap_vars(r, HAP2)->poss[0]);
    EXPECT_EQ(GT_REF_ALT1, hap_vars(r, HAP2)->orig_gts[0]);
    EXPECT_TRUE(logged(r, "1 overlapping variants in QUERY VCF, skipped"));
}

// The same downgrade applies in the other direction: a preceding deletion on haplotype 2 drops the
// haplotype-2 copy, leaving the haplotype-1 copy as 1|0.
TEST_F(ParseVariants, HomozygousDowngradedWhenHap2Overlaps) {
    ParseResult r = parse_records(dir, {record(100, "AGGT", "A", "0|1"),
                                        record(102, "C", "T", "1|1")});
    ASSERT_EQ(1, kept_on_hap(r, HAP1)); // the SNP only; the deletion is on the other haplotype
    ASSERT_EQ(1, kept_on_hap(r, HAP2)); // the deletion only; its SNP half overlapped
    EXPECT_EQ(101, hap_vars(r, HAP1)->poss[0]);
    EXPECT_EQ(GT_ALT1_REF, hap_vars(r, HAP1)->orig_gts[0]);
    EXPECT_TRUE(logged(r, "1 overlapping variants in QUERY VCF, skipped"));
}

// An insertion consumes no reference, so a colocated insertion on the other haplotype does not
// overlap it; the two-insertions rule drops it anyway, and the downgrade must follow that rule too.
TEST_F(ParseVariants, HomozygousInsertionDowngradedWhenColocatedWithInsertion) {
    ParseResult hap1_kept = parse_records(dir, {record(100, "A", "AGG", "0|1"),
                                                record(100, "A", "ATT", "1|1")});
    ASSERT_EQ(1, kept_on_hap(hap1_kept, HAP1));
    EXPECT_EQ("TT", hap_vars(hap1_kept, HAP1)->alts[0]);
    EXPECT_EQ(GT_ALT1_REF, hap_vars(hap1_kept, HAP1)->orig_gts[0]);
    EXPECT_TRUE(logged(hap1_kept, "1 overlapping variants in QUERY VCF, skipped"));

    ParseResult hap2_kept = parse_records(dir, {record(100, "A", "AGG", "1|0"),
                                                record(100, "A", "ATT", "1|1")});
    ASSERT_EQ(1, kept_on_hap(hap2_kept, HAP2));
    EXPECT_EQ("TT", hap_vars(hap2_kept, HAP2)->alts[0]);
    EXPECT_EQ(GT_REF_ALT1, hap_vars(hap2_kept, HAP2)->orig_gts[0]);
    EXPECT_TRUE(logged(hap2_kept, "1 overlapping variants in QUERY VCF, skipped"));
}

// Heterozygotes should land on either haplotype about equally; a lopsided split suggests the VCF
// was never really phased.
TEST_F(ParseVariants, HeterozygousImbalanceWarns) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "0|1"),
                                        record(200, "A", "G", "0|1"),
                                        record(300, "A", "G", "0|1")});
    EXPECT_TRUE(logged(r, "Imbalance of heterozygous variant phasing"));
}

TEST_F(ParseVariants, BalancedHeterozygotesDoNotWarn) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "0|1"),
                                        record(200, "A", "G", "1|0"),
                                        record(300, "A", "G", "0|1"),
                                        record(400, "A", "G", "1|0")});
    EXPECT_FALSE(logged(r, "Imbalance of heterozygous variant phasing"));
}

} // namespace
