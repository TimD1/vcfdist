/**
 * @file test_variant.cpp
 * @brief Unit tests for ctgVariants per-variant provenance and ploidy vectors.
 */
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include "gtest/gtest.h"
#include "../../../src/variant.h"
#include "../../../src/globals.h"

namespace {

/**************************************************************************************************/

// Single-sample header declaring the two contigs used by the fixtures below.
const std::string VCF_HEADER =
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=ctg1,length=200>\n"
    "##contig=<ID=ctgX,length=200>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n";

// Record ordinals within FIXTURE_BODY, pinned by the rec_idxs assertions below.
const int REC_HOM_SNP    = 0;
const int REC_HET_SNP    = 1;
const int REC_MULTIALLIC = 2;
const int REC_CPX        = 3;
const int REC_HAPLOID    = 4;

// One record per provenance case, in file order; ctgX carries the haploid call so that its
// ploidy of 1 does not conflict with the diploid ploidy recorded for ctg1.
const std::string FIXTURE_BODY =
    "ctg1\t11\t.\tA\tG\t30\tPASS\t.\tGT\t1|1\n"    // hom SNP, both haps, ALT 1
    "ctg1\t21\t.\tC\tT\t30\tPASS\t.\tGT\t0|1\n"    // het SNP, HAP2 only, ALT 1
    "ctg1\t31\t.\tA\tG,T\t30\tPASS\t.\tGT\t1|2\n"  // multiallelic, ALT 1 on HAP1, ALT 2 on HAP2
    "ctg1\t41\t.\tATTT\tG,GG\t30\tPASS\t.\tGT\t0|2\n" // CPX from ALT 2, split into INS+DEL
    "ctgX\t11\t.\tA\tC\t30\tPASS\t.\tGT\t1\n";    // haploid SNP, HAP1 only, ALT 1

/**
 * @brief Writes VCF text to a unique temporary file.
 * @param[in] contents Full VCF text, including header
 * @return Path to the newly written temporary file
 * @throws std::runtime_error Temporary file could not be created
 */
std::string write_temp_vcf(const std::string & contents) {
    std::string path_template =
        (std::filesystem::temp_directory_path() / "vcfdist_test_XXXXXX").string();
    std::vector<char> path(path_template.begin(), path_template.end());
    path.push_back('\0');
    int fd = mkstemp(path.data());
    if (fd < 0) throw std::runtime_error("failed to create temporary VCF");
    close(fd);
    std::ofstream out(path.data());
    out << contents;
    out.close();
    return std::string(path.data());
}

/**
 * @brief Parses a VCF fixture into a variantData container, then removes the temporary file.
 * @param[in] contents Full VCF text, including header
 * @return Populated variantData container for the QUERY callset
 */
std::shared_ptr<variantData> parse_fixture(const std::string & contents) {
    std::string vcf_fn = write_temp_vcf(contents);
    std::shared_ptr<variantData> vcf_data(new variantData());
    // parse_variants only stores the reference pointer, so no FASTA fixture is needed
    parse_variants(vcf_fn, vcf_data, std::shared_ptr<fastaData>(nullptr), QUERY);
    std::filesystem::remove(vcf_fn);
    return vcf_data;
}

/**************************************************************************************************/

/**
 * @class ProvenanceVectors
 * @brief Parses the shared VCF fixture once per test, with global state saved and restored.
 */
class ProvenanceVectors : public ::testing::Test {
protected:
    /** @brief Silences parser output and parses the shared fixture. */
    void SetUp() override {
        saved_verbosity = g.verbosity;
        g.verbosity = 0;
        vcf_data = parse_fixture(VCF_HEADER + FIXTURE_BODY);
        hap1 = vcf_data->variants[HAP1]["ctg1"];
        hap2 = vcf_data->variants[HAP2]["ctg1"];
        hapx = vcf_data->variants[HAP1]["ctgX"];
    }

    /** @brief Restores the global verbosity modified by SetUp. */
    void TearDown() override { g.verbosity = saved_verbosity; }

    int saved_verbosity = 1;
    std::shared_ptr<variantData> vcf_data;
    std::shared_ptr<ctgVariants> hap1;
    std::shared_ptr<ctgVariants> hap2;
    std::shared_ptr<ctgVariants> hapx;
};

/**************************************************************************************************/

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
    int ins = 3;
    int del = 4;

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
    std::shared_ptr<ctgVariants> vars(new ctgVariants("ctg1"));
    vars->add_var(10, 1, TYPE_SUB, BED_INSIDE, "A", "G", GT_ALT1_ALT1, 60, 60, 0);

    ASSERT_EQ(1, vars->n);
    EXPECT_EQ(-1, vars->rec_idxs[0]);
    EXPECT_EQ(-1, vars->alt_idxs[0]);
    EXPECT_EQ(0, vars->ploidies[0]);
}

/**************************************************************************************************/

} // namespace
