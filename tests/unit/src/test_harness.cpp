/**
 * @file test_harness.cpp
 * @brief Unit tests for test_helpers.cpp: the unit-test scaffolding itself.
 */
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/defs.h"
#include "../../../src/fasta.h"
#include "../../../src/globals.h"
#include "../../../src/variant.h"
#include "test_helpers.h"

namespace {

/* TempDir ****************************************************************************************/

TEST(TempDir, CreatesWritable) {
    GlobalsGuard guard;
    TempDir dir;
    ASSERT_TRUE(std::filesystem::is_directory(dir.path()));

    std::ofstream out(dir.path("probe.txt"));
    ASSERT_TRUE(out.is_open());
    out << "probe\n";
    out.close();
    EXPECT_TRUE(std::filesystem::is_regular_file(dir.path("probe.txt")));
}

TEST(TempDir, DistinctPaths) {
    GlobalsGuard guard;
    TempDir first;
    TempDir second;

    // mkdtemp() exists to guarantee this: two live instances never collide
    EXPECT_NE(first.path(), second.path());
    EXPECT_TRUE(std::filesystem::is_directory(first.path()));
    EXPECT_TRUE(std::filesystem::is_directory(second.path()));
}

TEST(TempDir, RemovesRecursively) {
    GlobalsGuard guard;
    std::string saved_path;
    {
        TempDir dir;
        saved_path = dir.path();
        ASSERT_TRUE(std::filesystem::create_directory(dir.path("nested")));
        std::ofstream out(dir.path("nested/probe.txt"));
        ASSERT_TRUE(out.is_open());
        out << "probe\n";
        out.close();
        ASSERT_TRUE(std::filesystem::exists(dir.path("nested/probe.txt")));
    }
    EXPECT_FALSE(std::filesystem::exists(saved_path));
}

TEST(TempDir, PathJoin) {
    GlobalsGuard guard;
    TempDir dir;

    // path() has no trailing slash, and path(name) adds exactly one separator
    EXPECT_NE('/', dir.path().back());
    EXPECT_EQ(dir.path() + "/probe.txt", dir.path("probe.txt"));
}

/* GlobalsGuard ***********************************************************************************/

TEST(GlobalsGuard, RestoresScalar) {
    const int before = g.sv_threshold;
    {
        GlobalsGuard guard;
        g.sv_threshold = before + 1234;
        ASSERT_EQ(before + 1234, g.sv_threshold);
    }
    EXPECT_EQ(before, g.sv_threshold);
}

TEST(GlobalsGuard, RestoresContainer) {
    const std::vector<std::string> before = g.filters;
    {
        GlobalsGuard guard;
        g.filters = {"PASS", "GT_LOWQUAL"};
        g.filter_ids = {0, 1};
        ASSERT_EQ(size_t(2), g.filters.size());
    }
    EXPECT_EQ(before, g.filters);
}

TEST(GlobalsGuard, SilencesVerbosity) {
    const int before = g.verbosity;
    {
        GlobalsGuard guard;
        EXPECT_EQ(0, g.verbosity);
    }
    EXPECT_EQ(before, g.verbosity);
}

/* write_tmp_vcf **********************************************************************************/

TEST(WriteTmpVcf, ParsesUnderHtslib) {
    GlobalsGuard guard;
    TempDir dir;
    std::string vcf_fn = write_tmp_vcf(dir,
            {"chr1\t5\t.\tA\tG\t60\tPASS\t.\tGT:GQ:PS\t0|1:60:1",
             "chr1\t9\t.\tA\tAC\t60\tPASS\t.\tGT:GQ:PS\t1|1:60:1"});

    std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
    std::shared_ptr<variantData> vcf(new variantData());
    parse_variants(vcf_fn, vcf, ref, TRUTH);

    EXPECT_EQ("TRUTH", vcf->sample);
    ASSERT_EQ(size_t(1), vcf->contigs.size());
    EXPECT_EQ("chr1", vcf->contigs[0]);
    EXPECT_EQ(12, vcf->lengths[0]);

    // the het SNP lands on one haplotype and the hom insertion on both
    std::shared_ptr<ctgVariants> hap1 = vcf->variants[HAP1]["chr1"];
    std::shared_ptr<ctgVariants> hap2 = vcf->variants[HAP2]["chr1"];
    EXPECT_EQ(3, hap1->n + hap2->n);
    std::shared_ptr<ctgVariants> snp_hap = (hap1->n == 2) ? hap1 : hap2;
    ASSERT_EQ(2, snp_hap->n);
    EXPECT_EQ(4, snp_hap->poss[0]);
    EXPECT_EQ("A", snp_hap->refs[0]);
    EXPECT_EQ("G", snp_hap->alts[0]);
    EXPECT_EQ(TYPE_SUB, snp_hap->types[0]);
    EXPECT_EQ(TYPE_INS, snp_hap->types[1]);
    EXPECT_EQ("C", snp_hap->alts[1]);
}

/* make_fasta *************************************************************************************/

TEST(MakeFasta, SequenceReadable) {
    GlobalsGuard guard;
    std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
    ASSERT_EQ(size_t(1), ref->fasta.size());
    EXPECT_EQ("ACGTACGTACGT", ref->fasta.at("chr1"));
    EXPECT_EQ(12, ref->lengths.at("chr1"));

    // lowercase input is uppercased by default, as the FASTA-parsing constructor does
    std::shared_ptr<fastaData> upper = make_fasta("chr1", "acgt");
    EXPECT_EQ("ACGT", upper->fasta.at("chr1"));
    std::shared_ptr<fastaData> verbatim = make_fasta("chr1", "acgt", false);
    EXPECT_EQ("acgt", verbatim->fasta.at("chr1"));

    // multi-contig sequences are stored independently
    std::shared_ptr<fastaData> multi = make_fasta({{"chr1", "ACGT"}, {"chrX", "TTGCAA"}});
    ASSERT_EQ(size_t(2), multi->fasta.size());
    EXPECT_EQ("ACGT", multi->fasta.at("chr1"));
    EXPECT_EQ("TTGCAA", multi->fasta.at("chrX"));
    EXPECT_EQ(6, multi->lengths.at("chrX"));
}

/* make_ctgVariants *******************************************************************************/

TEST(MakeCtgVariants, Roundtrip) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> vars = make_ctgVariants("chr1", {
            {4, 1, TYPE_SUB, "A", "G", GT_REF_ALT1, 40, 7, -1, BED_INSIDE},
            {8, 0, TYPE_INS, "", "CC", GT_ALT1_ALT1, 60, 7, 2, BED_OUTSIDE}});

    EXPECT_EQ("chr1", vars->ctg);
    ASSERT_EQ(2, vars->n);

    EXPECT_EQ(4, vars->poss[0]);
    EXPECT_EQ(1, vars->rlens[0]);
    EXPECT_EQ(TYPE_SUB, vars->types[0]);
    EXPECT_EQ(BED_INSIDE, vars->locs[0]);
    EXPECT_EQ("A", vars->refs[0]);
    EXPECT_EQ("G", vars->alts[0]);
    EXPECT_EQ(GT_REF_ALT1, vars->orig_gts[0]);
    EXPECT_FLOAT_EQ(40, vars->var_quals[0]);
    EXPECT_FLOAT_EQ(40, vars->gt_quals[0]);
    EXPECT_EQ(7, vars->phase_sets[0]);
    EXPECT_EQ(-1, vars->superclusters[0]);

    // order is preserved, and per-variant fields are not shared between entries
    EXPECT_EQ(8, vars->poss[1]);
    EXPECT_EQ(TYPE_INS, vars->types[1]);
    EXPECT_EQ("CC", vars->alts[1]);
    EXPECT_EQ(BED_OUTSIDE, vars->locs[1]);
    EXPECT_EQ(2, vars->superclusters[1]);

    // the per-haplotype lanes are sized alongside the parsed data
    ASSERT_EQ(size_t(HAPS), vars->errtypes.size());
    EXPECT_EQ(size_t(2), vars->errtypes[HAP1].size());
    EXPECT_EQ(size_t(2), vars->credit[HAP2].size());
}

TEST(MakeCtgVariants, ProvenanceFieldsReachTheirOwnVectors) {
    GlobalsGuard guard;

    // add_var() defaults every parameter past phase_set, so an omitted argument shifts the
    // remainder along silently. Distinct values are what catch that: a shift of even one position
    // lands one of these in a neighbouring vector.
    var_desc var;
    var.pos = 3;
    var.rlen = 1;
    var.ref = "T";
    var.alt = "C";
    var.phase_set = 11;
    var.supercluster = 5;
    var.rec_idx = 9;
    var.alt_idx = 2;
    var.ploidy = 1;
    std::shared_ptr<ctgVariants> vars = make_ctgVariants("chr1", {var});

    ASSERT_EQ(1, vars->n);
    EXPECT_EQ(11, vars->phase_sets[0]);
    EXPECT_EQ(9, vars->rec_idxs[0]);
    EXPECT_EQ(2, vars->alt_idxs[0]);
    EXPECT_EQ(1, vars->ploidies[0]);
    EXPECT_EQ(5, vars->superclusters[0]);

    // an unspecified field yields the same "unknown" sentinel add_var() defaults to
    std::shared_ptr<ctgVariants> plain = make_ctgVariants("chr1", {{4, 1, TYPE_SUB, "A", "G"}});
    EXPECT_EQ(-1, plain->rec_idxs[0]);
    EXPECT_EQ(-1, plain->alt_idxs[0]);
    EXPECT_EQ(0, plain->ploidies[0]);
    EXPECT_EQ(-1, plain->superclusters[0]);
}

/* alloc_reach_offs *******************************************************************************/

TEST(AllocReachOffs, SizeAndInit) {
    GlobalsGuard guard;

    // wf_swg_max_reach() indexes MATS matrices by (score % (max(x,o+e)+1)) and diagonal
    const int qlen = 5, tlen = 3, x = 5, o = 6, e = 2;
    std::vector<int> offs = alloc_reach_offs(qlen, tlen, x, o, e);
    ASSERT_EQ(size_t(MATS * (std::max(x, o+e) + 1) * (qlen + tlen - 1)), offs.size());
    for (size_t i = 0; i < offs.size(); i++) {
        ASSERT_EQ(-2, offs[i]) << "offset " << i;
    }

    // the score modulus is max(x, o+e)+1, so a large substitution penalty widens the buffer
    std::vector<int> sub_dominant = alloc_reach_offs(qlen, tlen, 20, o, e);
    EXPECT_EQ(size_t(MATS * 21 * (qlen + tlen - 1)), sub_dominant.size());
}

} // namespace
