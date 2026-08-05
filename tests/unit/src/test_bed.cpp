/**
 * @file test_bed.cpp
 * @brief Unit tests for bed.cpp: BED parsing, validation, and interval classification.
 */
#include <memory>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/bed.h"
#include "../../../src/defs.h"
#include "../../../src/fasta.h"
#include "../../../src/globals.h"
#include "../../../src/variant.h"
#include "test_helpers.h"

namespace {

/* Test helpers ***********************************************************************************/

/**
 * @brief Builds the two-region reference BED used by most bedData::contains tests.
 *
 * chr1 carries [10, 20) and [30, 40), leaving a [20, 30) gap between them and unbounded regions
 * on either side, so that every branch of contains() is reachable with in-range coordinates.
 * @return BED with two regions on chr1
 */
bedData two_region_bed() {
    return make_bed("chr1", {{10, 20}, {30, 40}});
}

/**
 * @brief Builds the two-contig reference shared by the intersect_contigs tests.
 * @return Reference holding 12-base chr1 and chr2 sequences
 */
std::shared_ptr<fastaData> two_contig_ref() {
    return make_fasta({{"chr1", "ACGTACGTACGT"}, {"chr2", "TTTTTTTTTTTT"}});
}

/* bedData::bedData *******************************************************************************/

TEST(BedCtor, OpenFailErrors) {
    GlobalsGuard guard;
    TempDir dir;
    EXPECT_EXIT(bedData bed(dir.path("absent.bed")), testing::ExitedWithCode(1),
            "Failed to open BED file");
}

TEST(BedCtor, ParsesThreeCols) {
    GlobalsGuard guard;

    // tiny.bed holds exactly "chr1\t2\t8"
    bedData bed(data_path("tiny.bed"));

    ASSERT_EQ(size_t(1), bed.contigs.size());
    EXPECT_EQ("chr1", bed.contigs[0]);
    ASSERT_EQ(size_t(1), bed.regions.count("chr1"));
    ASSERT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(2, bed.regions["chr1"].starts[0]);
    EXPECT_EQ(8, bed.regions["chr1"].stops[0]);

    // BED intervals are half-open, so [2, 8) covers six bases
    EXPECT_EQ(6L, bed.size);
}

TEST(BedCtor, IgnoresExtraCols) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8\tname\t0\t+"});

    bedData bed(bed_fn);

    ASSERT_EQ(size_t(1), bed.contigs.size());
    EXPECT_EQ("chr1", bed.contigs[0]);
    ASSERT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(2, bed.regions["chr1"].starts[0]);
    EXPECT_EQ(8, bed.regions["chr1"].stops[0]);
    EXPECT_EQ(6L, bed.size);
}

TEST(BedCtor, RunsCheck) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8", "chr1\t5\t10"});

    // the constructor validates after loading, so a malformed file is fatal at construction
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1), "BED overlap detected");
}

/* bedData::add ***********************************************************************************/

TEST(BedAdd, NewContig) {
    GlobalsGuard guard;
    bedData bed;

    bed.add("chr1", 2, 8);

    ASSERT_EQ(size_t(1), bed.contigs.size());
    EXPECT_EQ("chr1", bed.contigs[0]);
    ASSERT_EQ(size_t(1), bed.regions.size());
    ASSERT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(2, bed.regions["chr1"].starts[0]);
    EXPECT_EQ(8, bed.regions["chr1"].stops[0]);
}

TEST(BedAdd, AppendsSameContig) {
    GlobalsGuard guard;
    bedData bed;

    bed.add("chr1", 2, 8);
    bed.add("chr1", 10, 14);

    // a second region on a known contig appends rather than re-registering the contig
    ASSERT_EQ(size_t(1), bed.contigs.size());
    ASSERT_EQ(size_t(1), bed.regions.size());
    ASSERT_EQ(2, bed.regions["chr1"].n);
    EXPECT_EQ(std::vector<int>({2, 10}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({8, 14}), bed.regions["chr1"].stops);
}

TEST(BedAdd, AccumulatesSize) {
    GlobalsGuard guard;
    bedData bed;
    EXPECT_EQ(0L, bed.size);

    bed.add("chr1", 2, 8);   // 6 bases
    bed.add("chrX", 10, 14); // 4 bases

    EXPECT_EQ(10L, bed.size);
    EXPECT_EQ(std::vector<std::string>({"chr1", "chrX"}), bed.contigs);
}

/* bedData::check *********************************************************************************/

TEST(BedCheck, FlippedErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                bedData bed = make_bed("chr1", {{8, 2}});
                bed.check();
            }, testing::ExitedWithCode(1), "BED region chr1:8-2 stop precedes start");
}

TEST(BedCheck, ZeroLengthErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                bedData bed = make_bed("chr1", {{5, 5}});
                bed.check();
            }, testing::ExitedWithCode(1), "BED region chr1:5-5 length zero");
}

TEST(BedCheck, UnsortedErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                bedData bed = make_bed("chr1", {{10, 20}, {1, 5}});
                bed.check();
            }, testing::ExitedWithCode(1), "BED is unsorted");
}

TEST(BedCheck, OverlapErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                bedData bed = make_bed("chr1", {{2, 8}, {5, 10}});
                bed.check();
            }, testing::ExitedWithCode(1), "BED overlap detected");
}

TEST(BedCheck, SharedBoundaryWarns) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{2, 5}, {5, 8}});

    // abutting regions are mergeable but not invalid, so check() returns normally
    testing::internal::CaptureStderr();
    bed.check();
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("[WARN")) << out;
    EXPECT_NE(std::string::npos, out.find("should be merged")) << out;
}

TEST(BedCheck, ValidPasses) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{2, 5}, {7, 9}});

    testing::internal::CaptureStderr();
    bed.check();
    EXPECT_EQ("", testing::internal::GetCapturedStderr());
}

/* bedData::contains ******************************************************************************/

TEST(BedContains, NoBedInside) {
    GlobalsGuard guard;
    g.bed_exists = false;
    bedData bed; // deliberately empty: the early return precedes every lookup

    EXPECT_EQ(BED_INSIDE, bed.contains("chrZ", 100, 200, TYPE_SUB));
}

TEST(BedContains, FlippedErrors) {
    GlobalsGuard guard;
    g.bed_exists = true;
    EXPECT_EXIT({
                bedData bed = two_region_bed();
                bed.contains("chr1", 8, 2, TYPE_SUB);
            }, testing::ExitedWithCode(1), "Invalid region chr1:8-2 in BED contains");
}

TEST(BedContains, UnknownContigOffctg) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_OFFCTG, bed.contains("chrZ", 12, 15, TYPE_SUB));
}

TEST(BedContains, BeforeAllOutside) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 2, 5, TYPE_SUB));

    // the first region starts at 10, and a variant ending exactly there is still outside it
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 5, 10, TYPE_SUB));
}

TEST(BedContains, AfterAllOutside) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 45, 50, TYPE_SUB));

    // the last region ends at 40, and a variant starting exactly there is already past it
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 40, 45, TYPE_SUB));
}

TEST(BedContains, MiddleInside) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 12, 15, TYPE_SUB));
    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 32, 35, TYPE_SUB));

    // a variant filling a region exactly is inside it
    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 10, 20, TYPE_SUB));
}

TEST(BedContains, InsAtRegionEndBorder) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = two_region_bed();

    // an insertion anchored on the last base of [10, 20) adds sequence at the region edge, so it
    // is treated as straddling the boundary rather than contained
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 19, 19, TYPE_INS));

    // the special case is type-specific and edge-specific
    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 19, 19, TYPE_SUB));
    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 18, 18, TYPE_INS));
}

TEST(BedContains, BetweenRegionsOutside) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = two_region_bed();

    // wholly within the [20, 30) gap
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 22, 25, TYPE_SUB));

    // flush against both gap edges, which are exclusive on the left and inclusive on the right
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 20, 30, TYPE_SUB));
}

TEST(BedContains, PartialOverlapBorder) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = two_region_bed();

    // overhangs the left edge of the first region, so the start index falls off the front
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 5, 15, TYPE_SUB));

    // overhangs the right edge of the last region, so the stop index falls off the end
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 35, 45, TYPE_SUB));

    // straddles the right edge of the first region into the gap
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 18, 22, TYPE_SUB));

    // straddles the left edge of the second region out of the gap
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 28, 32, TYPE_SUB));
}

TEST(BedContains, SpansMultipleBorder) {
    GlobalsGuard guard;
    g.bed_exists = true;
    bedData bed = make_bed("chr1", {{10, 20}, {30, 40}, {50, 60}});

    // covers all of the middle region plus part of the outer two
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 15, 55, TYPE_SUB));
}

/* bedData::operator std::string ******************************************************************/

TEST(BedToString, RendersRegionsByContig) {
    GlobalsGuard guard;
    bedData bed;
    bed.add("chr1", 2, 8);
    bed.add("chr1", 10, 14);
    bed.add("chrX", 0, 5);

    // regions are stored in an unordered_map, so only per-contig content is deterministic
    std::string out = std::string(bed);
    EXPECT_NE(std::string::npos, out.find("chr1:\n\t2-8\n\t10-14\n")) << out;
    EXPECT_NE(std::string::npos, out.find("chrX:\n\t0-5\n")) << out;

    EXPECT_EQ("", std::string(bedData()));
}

/* intersect_contigs ******************************************************************************/

TEST(IntersectContigs, BedDropsExtraneous) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{0, 10}});
    std::shared_ptr<variantData> query =
            make_variantData(QUERY, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
    std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    testing::internal::GetCapturedStderr();

    // chr2 is absent from the BED, so it is dropped from the query along with its parallel fields
    EXPECT_EQ(std::vector<std::string>({"chr1"}), query->contigs);
    EXPECT_EQ(std::vector<int>({12}), query->lengths);
    EXPECT_EQ(std::vector< std::set<int> >({{2}}), query->observed_ploidies);
    EXPECT_EQ(size_t(0), query->variants[HAP1].count("chr2"));
    EXPECT_EQ(size_t(0), query->variants[HAP2].count("chr2"));

    // the reference is pruned to the BED contigs too
    EXPECT_EQ(size_t(1), ref->fasta.size());
    EXPECT_EQ(size_t(1), ref->fasta.count("chr1"));
}

TEST(IntersectContigs, BedMissingInFastaErrors) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed({{"chr1", {{0, 10}}}, {"chr2", {{0, 10}}}});
    EXPECT_EXIT({
                std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
                std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
                std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
                intersect_contigs(query, truth, ref);
            }, testing::ExitedWithCode(1), "Contig 'chr2' found in BED but not reference FASTA");
}

TEST(IntersectContigs, BedAddsEmptyContig) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed({{"chr1", {{0, 10}}}, {"chr2", {{0, 10}}}});
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
    std::shared_ptr<variantData> truth =
            make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    // a BED contig missing from the query is injected empty, having observed no ploidy at all
    EXPECT_EQ(std::vector<std::string>({"chr1", "chr2"}), query->contigs);
    EXPECT_EQ(std::vector<int>({12, 12}), query->lengths);
    EXPECT_EQ(std::vector< std::set<int> >({{2}, {}}), query->observed_ploidies);
    ASSERT_EQ(size_t(1), query->variants[HAP1].count("chr2"));
    ASSERT_EQ(size_t(1), query->variants[HAP2].count("chr2"));
    EXPECT_EQ(0, query->variants[HAP1]["chr2"]->n);
    EXPECT_EQ(0, query->variants[HAP2]["chr2"]->n);
    EXPECT_EQ("chr2", query->variants[HAP1]["chr2"]->ctg);

    EXPECT_NE(std::string::npos, out.find("found in BED but not query VCF")) << out;
}

// A contig injected empty has no observed ploidy at all, which is not a disagreement with the
// other callset. Matches "has ploidy" and "has ploidies" alike, so it pins both message forms.
TEST(IntersectContigs, InjectedEmptyContigDoesNotWarnOnPloidy) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed({{"chr1", {{0, 10}}}, {"chr2", {{0, 10}}}});
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
    std::shared_ptr<variantData> truth =
            make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_EQ(std::string::npos, out.find("contig 'chr2' has ploid")) << out;
}

TEST(IntersectContigs, NobedQueryOnlyContigWarns) {
    GlobalsGuard guard;
    g.bed_exists = false;
    std::shared_ptr<variantData> query =
            make_variantData(QUERY, {"chr1", "chr2"}, {12, 12}, {{2}, {1}});
    std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("[WARN")) << out;
    EXPECT_NE(std::string::npos,
            out.find("Contig 'chr2' found in query VCF but not truth VCF")) << out;

    // the truth gains an empty chr2 with no observed ploidy, so all query calls there are FPs
    EXPECT_EQ(std::vector<std::string>({"chr1", "chr2"}), truth->contigs);
    EXPECT_EQ(std::vector<int>({12, 12}), truth->lengths);
    EXPECT_EQ(std::vector< std::set<int> >({{2}, {}}), truth->observed_ploidies);
    ASSERT_EQ(size_t(1), truth->variants[HAP1].count("chr2"));
    EXPECT_EQ(0, truth->variants[HAP1]["chr2"]->n);
    EXPECT_EQ(0, truth->variants[HAP2]["chr2"]->n);
}

TEST(IntersectContigs, NobedTruthOnlyContigWarns) {
    GlobalsGuard guard;
    g.bed_exists = false;
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
    std::shared_ptr<variantData> truth =
            make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {1}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("[WARN")) << out;
    EXPECT_NE(std::string::npos,
            out.find("Contig 'chr2' found in truth VCF but not query VCF")) << out;

    // the query gains an empty chr2 with no observed ploidy, so all truth calls there are FNs
    EXPECT_EQ(std::vector<std::string>({"chr1", "chr2"}), query->contigs);
    EXPECT_EQ(std::vector<int>({12, 12}), query->lengths);
    EXPECT_EQ(std::vector< std::set<int> >({{2}, {}}), query->observed_ploidies);
    ASSERT_EQ(size_t(1), query->variants[HAP1].count("chr2"));
    EXPECT_EQ(0, query->variants[HAP1]["chr2"]->n);
    EXPECT_EQ(0, query->variants[HAP2]["chr2"]->n);
}

TEST(IntersectContigs, NobedFastaMissingErrors) {
    GlobalsGuard guard;
    g.bed_exists = false;
    EXPECT_EXIT({
                std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
                std::shared_ptr<variantData> truth =
                        make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
                std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
                intersect_contigs(query, truth, ref);
            }, testing::ExitedWithCode(1),
            "Contig 'chr2' found in truth VCF but not reference FASTA");
}

TEST(IntersectContigs, PloidyMismatchWarns) {
    GlobalsGuard guard;
    g.bed_exists = false;
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{1}});
    std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
    std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("[WARN")) << out;
    EXPECT_NE(std::string::npos,
            out.find("TRUTH contig 'chr1' has ploidies {2} and"
                     " QUERY contig 'chr1' has ploidies {1}"))
            << out;
}

// A contig legitimately carrying both ploidies in both callsets is not a disagreement.
TEST(IntersectContigs, MatchingMixedPloidyDoesNotWarn) {
    GlobalsGuard guard;
    g.bed_exists = false;
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{1, 2}});
    std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{1, 2}});
    std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_EQ(std::string::npos, out.find("has ploidies")) << out;
}

// Sets that overlap but are not equal still disagree: one callset saw a ploidy the other never did.
TEST(IntersectContigs, PartiallyOverlappingPloidySetsWarn) {
    GlobalsGuard guard;
    g.bed_exists = false;
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
    std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{1, 2}});
    std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos,
            out.find("TRUTH contig 'chr1' has ploidies {1,2} and"
                     " QUERY contig 'chr1' has ploidies {2}"))
            << out;
}

} // namespace
