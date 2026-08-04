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

/* Local helpers **********************************************************************************/

/**
 * @brief Builds a one-variant container with the given original and calculated genotypes.
 * @param[in] orig_gt Original genotype (GT_*) stored in orig_gts[0]
 * @param[in] calc_gt Calculated genotype (GT_*) stored in calc_gts[0]
 * @return Container holding a single SNP at chr1:100 with the requested genotypes
 */
std::shared_ptr<ctgVariants> make_gt_var(uint8_t orig_gt, uint8_t calc_gt) {
    var_desc var;
    var.pos = 100;
    var.rlen = 1;
    var.type = TYPE_SUB;
    var.ref = "A";
    var.alt = "C";
    var.gt = orig_gt;
    std::shared_ptr<ctgVariants> vars = make_ctgVariants("chr1", {var});
    vars->calc_gts[0] = calc_gt;
    return vars;
}

/**
 * @brief Builds a one-variant container of the given type with the given allele sequences.
 * @param[in] type Variant type (TYPE_*)
 * @param[in] ref Reference allele sequence
 * @param[in] alt Alternate allele sequence
 * @return Container holding a single variant at chr1:100
 */
std::shared_ptr<ctgVariants> make_typed_var(uint8_t type, const std::string & ref,
        const std::string & alt) {
    var_desc var;
    var.pos = 100;
    var.rlen = int(ref.size());
    var.type = type;
    var.ref = ref;
    var.alt = alt;
    return make_ctgVariants("chr1", {var});
}

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
    EXPECT_TRUE(vars.calc_gts.empty());
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
    ASSERT_EQ(size_t(PHASES), vars.errtypes.size());
    ASSERT_EQ(size_t(PHASES), vars.sync_group.size());
    ASSERT_EQ(size_t(PHASES), vars.callq.size());
    ASSERT_EQ(size_t(PHASES), vars.credit.size());
    ASSERT_EQ(size_t(PHASES), vars.ref_ed.size());
    ASSERT_EQ(size_t(PHASES), vars.query_ed.size());

    for (int hap = 0; hap < PHASES; hap++) {
        EXPECT_TRUE(vars.errtypes[hap].empty()) << "hap " << hap;
        EXPECT_TRUE(vars.sync_group[hap].empty()) << "hap " << hap;
        EXPECT_TRUE(vars.callq[hap].empty()) << "hap " << hap;
        EXPECT_TRUE(vars.credit[hap].empty()) << "hap " << hap;
        EXPECT_TRUE(vars.ref_ed[hap].empty()) << "hap " << hap;
        EXPECT_TRUE(vars.query_ed[hap].empty()) << "hap " << hap;
    }
}

TEST(CtgVariantsCtor, EmptyCtg) {
    GlobalsGuard guard;
    ctgVariants vars("");
    EXPECT_EQ("", vars.ctg);
    EXPECT_EQ(0, vars.n);
    EXPECT_EQ(size_t(PHASES), vars.errtypes.size());
}

/* add_var (copy overload) ************************************************************************/

// Asserts every per-variant field that exists today, so it must gain assertions as new per-variant
// vectors land, or "every field" quietly stops being every field: strata_bits (#47), is_phased
// (#46). rec_idxs/alt_idxs/ploidies (#48) are covered below.
TEST(AddVarCopy, Roundtrip) {
    GlobalsGuard guard;
    g.max_qual = 100; // above every quality used here, so no clamping obscures the copy
    std::shared_ptr<ctgVariants> src(new ctgVariants("chr20"));
    src->add_var(1234, 3, TYPE_DEL, BED_BORDER, "ACG", "A", GT_ALT1_REF, 44, 55, 77,
            101, 2, 1, 9, GT_REF_ALT1, ERRTYPE_TP, ERRTYPE_FP, 11, 12, 13.5, 14.5, 15, 16, 17, 18,
            0.25, 0.75);

    std::shared_ptr<ctgVariants> dst(new ctgVariants("chr20"));
    dst->add_var(src, 0);

    ASSERT_EQ(1, dst->n);
    EXPECT_EQ(1234, dst->poss[0]);
    EXPECT_EQ(3, dst->rlens[0]);
    EXPECT_EQ(TYPE_DEL, dst->types[0]);
    EXPECT_EQ(BED_BORDER, dst->locs[0]);
    EXPECT_EQ("ACG", dst->refs[0]);
    EXPECT_EQ("A", dst->alts[0]);
    EXPECT_EQ(GT_ALT1_REF, dst->orig_gts[0]);
    EXPECT_FLOAT_EQ(44, dst->gt_quals[0]);
    EXPECT_FLOAT_EQ(55, dst->var_quals[0]);
    EXPECT_EQ(77, dst->phase_sets[0]);
    EXPECT_EQ(101, dst->rec_idxs[0]);
    EXPECT_EQ(2, dst->alt_idxs[0]);
    EXPECT_EQ(1, dst->ploidies[0]);
    EXPECT_EQ(9, dst->superclusters[0]);
    EXPECT_EQ(GT_REF_ALT1, dst->calc_gts[0]);
    EXPECT_EQ(ERRTYPE_TP, dst->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_FP, dst->errtypes[HAP2][0]);
    EXPECT_EQ(11, dst->sync_group[HAP1][0]);
    EXPECT_EQ(12, dst->sync_group[HAP2][0]);
    EXPECT_FLOAT_EQ(13.5, dst->callq[HAP1][0]);
    EXPECT_FLOAT_EQ(14.5, dst->callq[HAP2][0]);
    EXPECT_EQ(15, dst->ref_ed[HAP1][0]);
    EXPECT_EQ(16, dst->ref_ed[HAP2][0]);
    EXPECT_EQ(17, dst->query_ed[HAP1][0]);
    EXPECT_EQ(18, dst->query_ed[HAP2][0]);
    EXPECT_FLOAT_EQ(0.25, dst->credit[HAP1][0]);
    EXPECT_FLOAT_EQ(0.75, dst->credit[HAP2][0]);

    // the copy overload forwards to add_var(), so the phasing lanes are re-defaulted, not copied
    EXPECT_EQ(PHASE_NONE, dst->phases[0]);
    EXPECT_EQ(PHASE_NONE, dst->pb_phases[0]);
    EXPECT_EQ(AC_UNKNOWN, dst->ac_errtype[0]);
}

TEST(AddVarCopy, PreservesHap1Hap2Distinct) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> src(new ctgVariants("chr20"));

    // every per-haplotype argument differs between haplotypes, so a transposition cannot pass
    src->add_var(100, 1, TYPE_SUB, BED_INSIDE, "A", "C", GT_ALT1_ALT1, 30, 30, 0,
            -1, -1, 2, 0, GT_ALT1_ALT1, ERRTYPE_TP, ERRTYPE_FN, 1, 2, 10, 20, 3, 4, 5, 6,
            0.1, 0.9);

    std::shared_ptr<ctgVariants> dst(new ctgVariants("chr20"));
    dst->add_var(src, 0);

    EXPECT_EQ(ERRTYPE_TP, dst->errtypes[HAP1][0]);
    EXPECT_EQ(ERRTYPE_FN, dst->errtypes[HAP2][0]);
    EXPECT_EQ(1, dst->sync_group[HAP1][0]);
    EXPECT_EQ(2, dst->sync_group[HAP2][0]);
    EXPECT_FLOAT_EQ(10, dst->callq[HAP1][0]);
    EXPECT_FLOAT_EQ(20, dst->callq[HAP2][0]);
    EXPECT_EQ(3, dst->ref_ed[HAP1][0]);
    EXPECT_EQ(4, dst->ref_ed[HAP2][0]);
    EXPECT_EQ(5, dst->query_ed[HAP1][0]);
    EXPECT_EQ(6, dst->query_ed[HAP2][0]);
    EXPECT_FLOAT_EQ(0.1, dst->credit[HAP1][0]);
    EXPECT_FLOAT_EQ(0.9, dst->credit[HAP2][0]);
}

TEST(AddVarCopy, SecondOfTwo) {
    GlobalsGuard guard;
    var_desc first;
    first.pos = 100;
    first.rlen = 1;
    first.ref = "A";
    first.alt = "C";
    var_desc second;
    second.pos = 200;
    second.rlen = 1;
    second.ref = "G";
    second.alt = "T";
    second.gt = GT_ALT1_ALT1;
    std::shared_ptr<ctgVariants> src = make_ctgVariants("chr20", {first, second});

    std::shared_ptr<ctgVariants> dst(new ctgVariants("chr20"));
    dst->add_var(src, 1);

    ASSERT_EQ(1, dst->n);
    EXPECT_EQ(200, dst->poss[0]);
    EXPECT_EQ("G", dst->refs[0]);
    EXPECT_EQ("T", dst->alts[0]);
    EXPECT_EQ(GT_ALT1_ALT1, dst->orig_gts[0]);
}

TEST(AddVarCopy, AppendsNotOverwrites) {
    GlobalsGuard guard;
    var_desc existing;
    existing.pos = 50;
    existing.rlen = 1;
    existing.ref = "A";
    existing.alt = "C";
    std::shared_ptr<ctgVariants> dst = make_ctgVariants("chr20", {existing});

    var_desc incoming;
    incoming.pos = 300;
    incoming.rlen = 1;
    incoming.ref = "G";
    incoming.alt = "T";
    std::shared_ptr<ctgVariants> src = make_ctgVariants("chr20", {incoming});

    dst->add_var(src, 0);

    ASSERT_EQ(2, dst->n);
    EXPECT_EQ(50, dst->poss[0]);
    EXPECT_EQ("A", dst->refs[0]);
    EXPECT_EQ(300, dst->poss[1]);
    EXPECT_EQ("G", dst->refs[1]);
}

/* add_var (full overload) ************************************************************************/

TEST(AddVar, AllFields) {
    GlobalsGuard guard;
    g.max_qual = 100;
    ctgVariants vars("chr20");
    vars.add_var(500, 2, TYPE_CPX, BED_OUTSIDE, "AC", "GT", GT_ALT1_ALT1, 21, 22, 33,
            12, 3, 2, 7, GT_ALT1_REF, ERRTYPE_FN, ERRTYPE_TP, 4, 5, 6.5, 7.5, 8, 9, 10, 11,
            0.4, 0.6);

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
    EXPECT_EQ(GT_ALT1_REF, vars.calc_gts[0]);
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
    vars.add_var(100, 1, TYPE_SUB, BED_INSIDE, "A", "C", GT_REF_ALT1, 99, 99, 0);

    EXPECT_FLOAT_EQ(60, vars.var_quals[0]);

    // only var_qual is clamped; gt_qual is stored verbatim
    EXPECT_FLOAT_EQ(99, vars.gt_quals[0]);
}

TEST(AddVar, QualBelowCap) {
    GlobalsGuard guard;
    g.max_qual = 60;
    ctgVariants vars("chr20");
    vars.add_var(100, 1, TYPE_SUB, BED_INSIDE, "A", "C", GT_REF_ALT1, 30, 30, 0);
    EXPECT_FLOAT_EQ(30, vars.var_quals[0]);
}

TEST(AddVar, QualNegative) {
    GlobalsGuard guard;
    g.max_qual = 60;
    ctgVariants vars("chr20");
    vars.add_var(100, 1, TYPE_SUB, BED_INSIDE, "A", "C", GT_REF_ALT1, -5, -5, 0);

    // std::min() only caps from above, so a negative quality is stored unchanged
    EXPECT_FLOAT_EQ(-5, vars.var_quals[0]);
}

TEST(AddVar, PushesPhaseDefaults) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(100, 1, TYPE_SUB, BED_INSIDE, "A", "C", GT_REF_ALT1, 30, 30, 0);

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
    vars.add_var(100, 1, TYPE_SUB, BED_INSIDE, "A", "C", GT_REF_ALT1, 30, 30, 0);

    // the rec_idx/alt_idx/ploidy defaults are asserted by ProvenanceVectorDefaults.UnknownSentinels
    EXPECT_EQ(-1, vars.superclusters[0]);
    EXPECT_EQ(GT_REF_REF, vars.calc_gts[0]);
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

TEST(AddVar, LaneLengthsTrackN) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    for (int i = 0; i < 3; i++) {
        vars.add_var(100*i, 1, TYPE_SUB, BED_INSIDE, "A", "C", GT_REF_ALT1, 30, 30, 0);
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
    EXPECT_EQ(n, vars.calc_gts.size());
    EXPECT_EQ(n, vars.phases.size());
    EXPECT_EQ(n, vars.pb_phases.size());
    EXPECT_EQ(n, vars.ac_errtype.size());
    for (int hap = 0; hap < PHASES; hap++) {
        EXPECT_EQ(n, vars.errtypes[hap].size()) << "hap " << hap;
        EXPECT_EQ(n, vars.sync_group[hap].size()) << "hap " << hap;
        EXPECT_EQ(n, vars.callq[hap].size()) << "hap " << hap;
        EXPECT_EQ(n, vars.ref_ed[hap].size()) << "hap " << hap;
        EXPECT_EQ(n, vars.query_ed[hap].size()) << "hap " << hap;
        EXPECT_EQ(n, vars.credit[hap].size()) << "hap " << hap;
    }
}

TEST(AddVar, InsRlenZero) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(100, 0, TYPE_INS, BED_INSIDE, "", "ACGT", GT_REF_ALT1, 30, 30, 0);

    ASSERT_EQ(1, vars.n);
    EXPECT_EQ(0, vars.rlens[0]);
    EXPECT_EQ(TYPE_INS, vars.types[0]);
    EXPECT_EQ("", vars.refs[0]);
    EXPECT_EQ("ACGT", vars.alts[0]);
}

TEST(AddVar, DelEmptyAlt) {
    GlobalsGuard guard;
    ctgVariants vars("chr20");
    vars.add_var(100, 4, TYPE_DEL, BED_INSIDE, "ACGT", "", GT_REF_ALT1, 30, 30, 0);

    ASSERT_EQ(1, vars.n);
    EXPECT_EQ(4, vars.rlens[0]);
    EXPECT_EQ(TYPE_DEL, vars.types[0]);
    EXPECT_EQ("ACGT", vars.refs[0]);
    EXPECT_EQ("", vars.alts[0]);
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
        vcf_opts opts;
        opts.sample = "QUERY";
        opts.contigs = {"##contig=<ID=ctg1,length=200>", "##contig=<ID=ctgX,length=200>"};
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

// The copying overload must carry provenance and ploidy through unchanged.
TEST_F(ProvenanceVectors, CopyingOverloadPreservesProvenance) {
    std::shared_ptr<ctgVariants> copy(new ctgVariants("ctg1"));
    for (int vi = 0; vi < hap2->n; vi++) copy->add_var(hap2, vi);

    ASSERT_EQ(hap2->n, copy->n);
    EXPECT_EQ(hap2->rec_idxs, copy->rec_idxs);
    EXPECT_EQ(hap2->alt_idxs, copy->alt_idxs);
    EXPECT_EQ(hap2->ploidies, copy->ploidies);
}

// Callers with no source record (e.g. CIGAR-derived variants) get the unknown sentinels.
TEST(ProvenanceVectorDefaults, UnknownSentinels) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars(new ctgVariants("ctg1"));
    vars->add_var(10, 1, TYPE_SUB, BED_INSIDE, "A", "G", GT_ALT1_ALT1, 60, 60, 0);

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

TEST(SetAlleleErrtype, ZeroToOne) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_REF);
    EXPECT_EQ(AC_ERR_0_TO_1, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_0_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, ZeroToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EQ(AC_ERR_0_TO_2, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_0_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, OneToZero) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_REF_ALT1);
    EXPECT_EQ(AC_ERR_1_TO_0, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_1_TO_0, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, OneToOne) {
    GlobalsGuard guard;

    // a heterozygous call on the opposite haplotype is still one allele called for one expected
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_REF, GT_REF_ALT1);
    EXPECT_EQ(AC_ERR_1_TO_1, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_1_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, OneToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    EXPECT_EQ(AC_ERR_1_TO_2, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_1_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TwoToZero) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_0, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_2_TO_0, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TwoToOne) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_1, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_2_TO_1, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, TwoToTwo) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    EXPECT_EQ(AC_ERR_2_TO_2, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_ERR_2_TO_2, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, RefrefRefrefUnknown) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_REF_REF);

    // seed a different value so that a return-without-store would leave it behind
    vars->ac_errtype[0] = AC_ERR_2_TO_2;

    EXPECT_EQ(AC_UNKNOWN, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_UNKNOWN, vars->ac_errtype[0]);
}

TEST(SetAlleleErrtype, HaploidUnknown) {
    GlobalsGuard guard;

    // a haploid calc_gt matches none of the diploid branches, so it falls through to AC_UNKNOWN
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1, GT_ALT1);
    vars->ac_errtype[0] = AC_ERR_1_TO_1;

    EXPECT_EQ(AC_UNKNOWN, vars->set_allele_errtype(0));
    EXPECT_EQ(AC_UNKNOWN, vars->ac_errtype[0]);
}

/* calcgt_is_swapped ******************************************************************************/

TEST(CalcgtIsSwapped, EqualFalse) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_ALT1);
    EXPECT_FALSE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, OrigHomFalse) {
    GlobalsGuard guard;

    // orig 1|1 reports data from both haplotypes, so haplotype order does not matter
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    EXPECT_FALSE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, OrigRefFalse) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_REF, GT_REF_ALT1);
    EXPECT_FALSE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, CalcRefFalse) {
    GlobalsGuard guard;

    // calc 0|0 reports no data, so haplotype order does not matter
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_REF_REF);
    EXPECT_FALSE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, HetOpposite0110True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_REF);
    EXPECT_TRUE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, HetOpposite1001True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_REF, GT_REF_ALT1);
    EXPECT_TRUE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, Orig01CalcHomCreditHap1True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.9;
    vars->credit[HAP2][0] = 0.1;

    // orig 0|1 expects the allele on HAP2, so the better HAP1 credit is reported by swapping
    EXPECT_TRUE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, Orig01CalcHomCreditHap2False) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.1;
    vars->credit[HAP2][0] = 0.9;
    EXPECT_FALSE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, Orig10CalcHomCreditHap2True) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_REF, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.1;
    vars->credit[HAP2][0] = 0.9;
    EXPECT_TRUE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, CreditTieFalse) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_REF_ALT1, GT_ALT1_ALT1);
    vars->credit[HAP1][0] = 0.5;
    vars->credit[HAP2][0] = 0.5;

    // the comparison is a strict >, so equal credits keep the original haplotype order
    EXPECT_FALSE(vars->calcgt_is_swapped(0));
}

TEST(CalcgtIsSwapped, UnexpectedErrors) {
    GlobalsGuard guard;

    // a haploid orig_gt reaches no branch, so the final else reports an unexpected pair
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1, GT_REF_ALT1);
    EXPECT_EXIT(vars->calcgt_is_swapped(0), testing::ExitedWithCode(1),
            "Unexpected orig/calc genotypes for variant");
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

TEST(VarOnHap, HapGt1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EXIT(vars->var_on_hap(0, 2), testing::ExitedWithCode(1), "Unexpected haplotype 2");
}

/* set_var_calcgt_on_hap **************************************************************************/

TEST(SetVarCalcgtOnHap, RefrefSetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    vars->set_var_calcgt_on_hap(0, HAP1, true);
    EXPECT_EQ(GT_ALT1_REF, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, RefrefSetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    vars->set_var_calcgt_on_hap(0, HAP2, true);
    EXPECT_EQ(GT_REF_ALT1, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, RefrefUnsetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP1, false), testing::ExitedWithCode(1),
            "Variant calc_gt already unset");
}

TEST(SetVarCalcgtOnHap, RefrefUnsetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP2, false), testing::ExitedWithCode(1),
            "Variant calc_gt already unset");
}

TEST(SetVarCalcgtOnHap, RefaltSetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    vars->set_var_calcgt_on_hap(0, HAP1, true);
    EXPECT_EQ(GT_ALT1_ALT1, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, RefaltSetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP2, true), testing::ExitedWithCode(1),
            "Variant calc_gt already set");
}

TEST(SetVarCalcgtOnHap, RefaltUnsetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP1, false), testing::ExitedWithCode(1),
            "Variant calc_gt already unset");
}

TEST(SetVarCalcgtOnHap, RefaltUnsetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_ALT1);
    vars->set_var_calcgt_on_hap(0, HAP2, false);
    EXPECT_EQ(GT_REF_REF, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltrefSetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Variant calc_gt already set");
}

TEST(SetVarCalcgtOnHap, AltrefSetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    vars->set_var_calcgt_on_hap(0, HAP2, true);
    EXPECT_EQ(GT_ALT1_ALT1, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltrefUnsetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    vars->set_var_calcgt_on_hap(0, HAP1, false);
    EXPECT_EQ(GT_REF_REF, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltrefUnsetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_REF);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP2, false), testing::ExitedWithCode(1),
            "Variant calc_gt already unset");
}

TEST(SetVarCalcgtOnHap, AltaltSetHap1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Variant calc_gt already set");
}

TEST(SetVarCalcgtOnHap, AltaltSetHap2Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP2, true), testing::ExitedWithCode(1),
            "Variant calc_gt already set");
}

TEST(SetVarCalcgtOnHap, AltaltUnsetHap1) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);

    // clearing HAP1 leaves the alternate allele on HAP2 only
    vars->set_var_calcgt_on_hap(0, HAP1, false);
    EXPECT_EQ(GT_REF_ALT1, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, AltaltUnsetHap2) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1_ALT1);
    vars->set_var_calcgt_on_hap(0, HAP2, false);
    EXPECT_EQ(GT_ALT1_REF, vars->calc_gts[0]);
}

TEST(SetVarCalcgtOnHap, MissingErrors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_MISSING);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Unexpected calc_gts value");
}

TEST(SetVarCalcgtOnHap, HaploidErrors) {
    GlobalsGuard guard;

    // a haploid calc_gt matches none of the four diploid states, so the default branch errors
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_ALT1);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP1, true), testing::ExitedWithCode(1),
            "Unexpected calc_gts value");
}

TEST(SetVarCalcgtOnHap, HapGt1Errors) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_REF_REF);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, 2, true), testing::ExitedWithCode(1),
            "Unexpected hap idx 2");
}

TEST(SetVarCalcgtOnHap, ErrorsSuppressedWithIgnore) {
    GlobalsGuard guard;

    // every invalid (state, hap, set) triple from the four diploid states
    struct transition { uint8_t calc_gt; int hap; bool set; };
    const std::vector<transition> invalid = {
        {GT_REF_REF,   HAP1, false}, {GT_REF_REF,   HAP2, false},
        {GT_REF_ALT1,  HAP2, true},  {GT_REF_ALT1,  HAP1, false},
        {GT_ALT1_REF,  HAP1, true},  {GT_ALT1_REF,  HAP2, false},
        {GT_ALT1_ALT1, HAP1, true},  {GT_ALT1_ALT1, HAP2, true},
    };

    for (const transition & t : invalid) {
        std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, t.calc_gt);
        vars->set_var_calcgt_on_hap(0, t.hap, t.set, true);
        EXPECT_EQ(t.calc_gt, vars->calc_gts[0])
                << "calc_gt " << int(t.calc_gt) << " hap " << t.hap << " set " << t.set;
    }
}

TEST(SetVarCalcgtOnHap, IgnoreDoesNotSuppressMissing) {
    GlobalsGuard guard;

    // ignore_errors guards only the four diploid states; the default branch always errors
    std::shared_ptr<ctgVariants> vars = make_gt_var(GT_ALT1_ALT1, GT_MISSING);
    EXPECT_EXIT(vars->set_var_calcgt_on_hap(0, HAP1, true, true), testing::ExitedWithCode(1),
            "Unexpected calc_gts value");
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
    for (int hap = 0; hap < HAPS; hap++) {
        EXPECT_TRUE(vcf.variants[hap].empty()) << "hap " << hap;
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
    EXPECT_TRUE(vcf.ploidy.empty());
}

/* parse-time filtering, counters, and summary warnings *******************************************/

/**
 * @brief Counts variants that survived parsing on one haplotype of chr1.
 * @param[in] r Result of parse_records()
 * @param[in] hap Haplotype index (HAP1 or HAP2)
 * @return Number of surviving variants on that haplotype
 */
int kept_on_hap(const ParseResult & r, int hap) {
    return r.vars->variants[hap]["chr1"]->n;
}

/**
 * @brief Counts variants that survived parsing across both haplotypes of chr1.
 * @param[in] r Result of parse_records()
 * @return Total number of surviving variants
 */
int total_kept(const ParseResult & r) {
    return kept_on_hap(r, HAP1) + kept_on_hap(r, HAP2);
}

/**
 * @brief Reports whether the written VCF contains a record at a given position.
 * @param[in] r Result of parse_records()
 * @param[in] pos 1-based VCF position
 * @return True if a chr1 data line at that position was written
 */
bool wrote_pos(const ParseResult & r, int pos) {
    return r.out_vcf.find("\nchr1\t" + std::to_string(pos) + "\t") != std::string::npos;
}

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

/* reasons that stay drops ************************************************************************/

// A spanning deletion allele carries no variation, so it is dropped and counted.
TEST_F(ParseVariants, SpanningDeletionDroppedAndCounted) {
    ParseResult r = parse_records(dir, {record(100, "A", "*", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(wrote_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants spanned by deletion in QUERY VCF, skipped"));
}

// An ALT identical to its REF carries no variation, whether one base long or several.
TEST_F(ParseVariants, RefCallDroppedAndCounted) {
    ParseResult r = parse_records(dir, {record(100, "A", "A", "1|0"),
                                        record(200, "AT", "AT", "1|0")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(wrote_pos(r, 100));
    EXPECT_FALSE(wrote_pos(r, 200));
    EXPECT_TRUE(logged(r, "2 reference variants in QUERY VCF, skipped"));
}

// A no-call has no known allele on either haplotype, so the whole record is dropped.
TEST_F(ParseVariants, NoCallDroppedAndCounted) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", ".|.")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(wrote_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants with no known alleles (.|.) in QUERY VCF, skipped"));
    EXPECT_TRUE(logged(r, ".|.: 1"));
}

// A no-call is one dropped record, not one dropped record per haplotype.
TEST_F(ParseVariants, NoCallCountedOncePerRecord) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", ".|.")});
    EXPECT_TRUE(logged(r, "1 variants with no known alleles"));
    EXPECT_FALSE(logged(r, "2 variants with no known alleles"));
}

/* half calls *************************************************************************************/

// A half call keeps its known allele, so it must not be tallied as a no-call.
TEST_F(ParseVariants, HalfCallCountedDistinctlyFromNoCall) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|.")});
    EXPECT_TRUE(logged(r, "X|.: 1"));
    EXPECT_FALSE(logged(r, ".|.:")); // histogram line is only printed for nonzero counts
    EXPECT_FALSE(logged(r, "no known alleles"));
}

// The summary must not claim a half call was skipped, because its known allele was evaluated.
TEST_F(ParseVariants, HalfCallNotReportedAsSkipped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|.")});
    EXPECT_EQ(1, total_kept(r));
    EXPECT_TRUE(wrote_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants with a half call (1|.) in QUERY VCF, known allele kept"));
    EXPECT_FALSE(logged(r, "skipped"));
}

// A missing allele on either haplotype leaves the known allele on the other.
TEST_F(ParseVariants, HalfCallKeptOnTheHaplotypeWithTheKnownAllele) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|."),
                                        record(200, "A", "G", ".|1")});
    EXPECT_EQ(1, kept_on_hap(r, HAP1));
    EXPECT_EQ(1, kept_on_hap(r, HAP2));
    EXPECT_TRUE(logged(r, "X|.: 2"));
    EXPECT_TRUE(logged(r, "2 variants with a half call"));
}

/* reasons deferred to separate work **************************************************************/

// Overlapping variants are still dropped at parse time, with their counter and warning intact.
TEST_F(ParseVariants, OverlappingVariantStillDropped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "1|0"),
                                        record(100, "A", "T", "1|0")});
    EXPECT_EQ(1, kept_on_hap(r, HAP1));
    EXPECT_EQ(0, kept_on_hap(r, HAP2));
    EXPECT_EQ(std::string::npos, r.out_vcf.find("\tT\t")); // second, overlapping ALT
    EXPECT_TRUE(logged(r, "1 overlapping variants in QUERY VCF, skipped"));
}

// Unphased heterozygous genotypes are still dropped at parse time, counter and warning intact.
TEST_F(ParseVariants, UnphasedHeterozygousGenotypeStillDropped) {
    ParseResult r = parse_records(dir, {record(100, "A", "G", "0/1")});
    EXPECT_EQ(0, total_kept(r));
    EXPECT_FALSE(wrote_pos(r, 100));
    EXPECT_TRUE(logged(r, "1 variants with unphased genotypes in QUERY VCF, skipped"));
}

} // namespace
