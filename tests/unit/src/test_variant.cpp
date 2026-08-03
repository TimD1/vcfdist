/**
 * @file test_variant.cpp
 * @brief Unit tests for variant.cpp: genotype, allele-count, and variant-type logic.
 */
#include <fcntl.h>
#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/globals.h"
#include "../../../src/variant.h"
#include "test_helpers.h"

namespace {

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

/* parse-time filtering, counters, and summary warnings *******************************************/

/** @brief Parsed variants plus everything parse_variants() reported to stderr. */
struct ParseResult {
    std::shared_ptr<variantData> vars; ///< Variants that survived parse-time filtering
    std::string log;                   ///< All INFO/WARN output from parse_variants()
    std::string out_vcf;               ///< VCF written from the surviving variants
};

/**
 * @brief Builds a single-sample VCF record line on chr1 with phase set 1.
 * @param[in] pos 1-based VCF position
 * @param[in] ref REF allele
 * @param[in] alt ALT allele
 * @param[in] gt GT field value (e.g. "1|0", "1|.", ".|.")
 * @return One tab-separated VCF data line, without a trailing newline
 */
std::string record(int pos, const std::string & ref, const std::string & alt,
        const std::string & gt) {
    return "chr1\t" + std::to_string(pos) + "\t.\t" + ref + "\t" + alt +
        "\t50\tPASS\t.\tGT:PS\t" + gt + ":1";
}

/**
 * @brief Reads an entire file into a string.
 * @param[in] fn Input filename
 * @return File contents, or an empty string if the file cannot be opened
 */
std::string read_text(const std::string & fn) {
    std::ifstream in(fn);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

/**
 * @class StderrToFile
 * @brief Redirects the C `stderr` stream to a file for the object's lifetime.
 *
 * WARN() and INFO() reach stderr through fprintf(), so capturing the summary means redirecting
 * the underlying file descriptor; swapping std::cerr's streambuf would not intercept it.
 * Restoring in the destructor keeps a failed assertion from leaving stderr pointing into the
 * temporary directory after TempDir has deleted it.
 */
class StderrToFile {
public:
    /** @brief Redirects stderr to fn, truncating any existing contents. */
    explicit StderrToFile(const std::string & fn)
            : saved_fd(dup(fileno(stderr))),
              file_fd(open(fn.data(), O_WRONLY | O_CREAT | O_TRUNC, 0644)) {
        std::fflush(stderr);
        dup2(file_fd, fileno(stderr));
    }

    /** @brief Flushes the redirected output and restores the original stderr. */
    ~StderrToFile() {
        std::fflush(stderr);
        dup2(saved_fd, fileno(stderr));
        close(saved_fd);
        close(file_fd);
    }

    StderrToFile(const StderrToFile &) = delete;
    StderrToFile & operator=(const StderrToFile &) = delete;

private:
    int saved_fd; ///< Duplicate of the original stderr descriptor
    int file_fd;  ///< Descriptor of the redirect target
};

/**
 * @brief Parses VCF records with parse_variants(), capturing its stderr and output VCF.
 * @param[in] dir Temporary directory owning the fixture and captured output
 * @param[in] records VCF data lines, without trailing newlines
 * @return Surviving variants, captured log output, and the VCF written from those variants
 * @note The written VCF stands in for summary.vcf: both are generated from the variants that
 *       survive parse-time filtering, so a variant absent here is absent from summary.vcf.
 */
ParseResult parse_records(const TempDir & dir, const std::vector<std::string> & records) {
    vcf_opts opts;
    opts.sample = "QUERY";
    opts.contigs = {"##contig=<ID=chr1,length=1000>"};
    const std::string vcf_fn = write_tmp_vcf(dir, records, opts);
    const std::string log_fn = dir.path("parse.log");
    const std::string out_fn = dir.path("out.vcf");

    ParseResult result;
    result.vars = std::make_shared<variantData>();
    std::shared_ptr<fastaData> ref = make_fasta("chr1", std::string(1000, 'A'));

    { // stderr is redirected for the parse alone, so the INFO/WARN summary can be asserted on
        StderrToFile redirect(log_fn);
        parse_variants(vcf_fn, result.vars, ref, QUERY);
    }

    result.log = read_text(log_fn);
    result.vars->write_vcf(out_fn);
    result.out_vcf = read_text(out_fn);
    return result;
}

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
 * @brief Reports whether the log contains a substring.
 * @param[in] r Result of parse_records()
 * @param[in] text Substring to search for
 * @return True if the log contains the substring
 */
bool logged(const ParseResult & r, const std::string & text) {
    return r.log.find(text) != std::string::npos;
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
