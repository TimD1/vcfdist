/**
 * @file test_bed.cpp
 * @brief Unit tests for bed.cpp: BED parsing, validation, and interval classification.
 */
#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
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

/// Lines the three encoding tests share, so any difference between them is the encoding alone.
const std::vector<std::string> encoded_lines = {"chr1\t2\t8", "chr1\t20\t30", "chr2\t0\t5"};

/**
 * @brief Asserts that a bedData holds exactly the regions encoded_lines describes.
 * @param[in] bed BED parsed from encoded_lines in one of the three accepted encodings
 */
void expect_encoded_regions(bedData & bed) {
    ASSERT_EQ(size_t(2), bed.contigs.size());
    EXPECT_EQ("chr1", bed.contigs[0]);
    EXPECT_EQ("chr2", bed.contigs[1]);

    ASSERT_EQ(2, bed.regions["chr1"].n);
    EXPECT_EQ(2, bed.regions["chr1"].starts[0]);
    EXPECT_EQ(8, bed.regions["chr1"].stops[0]);
    EXPECT_EQ(20, bed.regions["chr1"].starts[1]);
    EXPECT_EQ(30, bed.regions["chr1"].stops[1]);

    ASSERT_EQ(1, bed.regions["chr2"].n);
    EXPECT_EQ(0, bed.regions["chr2"].starts[0]);
    EXPECT_EQ(5, bed.regions["chr2"].stops[0]);

    EXPECT_EQ(21L, bed.size);
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

// The same file BedCtor.RunsCheck rejects: normalizing repairs the overlap, so the check that
// still runs afterwards has nothing left to reject.
TEST(BedCtor, NormalizedOverlapPassesCheck) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8", "chr1\t5\t10"});

    testing::internal::CaptureStderr();
    bedData bed(bed_fn, true);
    EXPECT_EQ("", testing::internal::GetCapturedStderr());

    ASSERT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(2, bed.regions["chr1"].starts[0]);
    EXPECT_EQ(10, bed.regions["chr1"].stops[0]);
    EXPECT_EQ(8L, bed.size);
}

// Normalizing repairs disorder and overlap and nothing else, so the validation it precedes is not
// weakened: a malformation sorting and merging cannot fix is still fatal.
TEST(BedCtor, NormalizeStillRunsCheck) {
    GlobalsGuard guard;
    TempDir dir;

    EXPECT_EXIT({
                std::string bed_fn = write_tmp_bed(dir, {"chr1\t8\t2"}, "flipped.bed");
                bedData bed(bed_fn, true);
            }, testing::ExitedWithCode(1), "BED region chr1:8-2 stop precedes start");

    EXPECT_EXIT({
                std::string bed_fn = write_tmp_bed(dir, {"chr1\t5\t5"}, "empty_region.bed");
                bedData bed(bed_fn, true);
            }, testing::ExitedWithCode(1), "BED region chr1:5-5 length zero");
}

TEST(BedCtor, RecordsFilename) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8"});

    // the messages normalize() prints name the file, so the file has to be remembered
    bedData bed(bed_fn);

    EXPECT_EQ(bed_fn, bed.filename);
    EXPECT_EQ("", bedData().filename);
}

TEST(BedCtor, NonNumericCoordErrors) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8", "chr1\tstart\tstop"});

    // the offending field and its line are named, since the filename itself is fine here
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1),
            "Invalid coordinate 'start' on line 2 of BED file");
}

TEST(BedCtor, PartlyNumericCoordErrors) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8bp"});

    // a coordinate must be numeric in full; a trailing suffix is not silently dropped
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1),
            "Invalid coordinate '8bp' on line 1 of BED file");
}

TEST(BedCtor, OutOfRangeCoordErrors) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t99999999999"});

    // a coordinate too large for int is reported the same way as a non-numeric one
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1),
            "Invalid coordinate '99999999999' on line 1 of BED file");
}

TEST(BedCtor, MissingColumnErrors) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2"});

    // an absent stop column reads as an empty coordinate rather than parsing as zero
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1),
            "Invalid coordinate '' on line 1 of BED file");
}

TEST(BedCtor, BlankLineErrors) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8", ""});

    // a trailing blank line has no coordinates to read, so it is rejected with its line number
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1),
            "Invalid coordinate '' on line 2 of BED file");
}

TEST(BedCtor, ParsesUncompressed) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, encoded_lines, "test.bed", BEDZIP_NONE);

    bedData bed(bed_fn);

    expect_encoded_regions(bed);
}

TEST(BedCtor, ParsesGzip) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, encoded_lines, "test.bed.gz", BEDZIP_GZIP);

    // the GIAB stratification sets ship as plain gzip, which is not seekable and not blocked
    bedData bed(bed_fn);

    expect_encoded_regions(bed);
}

TEST(BedCtor, ParsesBgzip) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, encoded_lines, "test.bed.gz", BEDZIP_BGZIP);

    bedData bed(bed_fn);

    expect_encoded_regions(bed);
}

TEST(BedCtor, DetectsEncodingFromContents) {
    GlobalsGuard guard;
    TempDir dir;
    // the extension says plain text and the bytes say bgzip; the bytes are what must decide
    std::string bed_fn = write_tmp_bed(dir, encoded_lines, "misnamed.bed", BEDZIP_BGZIP);

    bedData bed(bed_fn);

    expect_encoded_regions(bed);
}

TEST(BedCtor, MalformedCoordInGzipNamesLine) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8", "chr1\tstart\tstop"},
            "test.bed.gz", BEDZIP_GZIP);

    // line numbers count decoded lines, so a compressed file reports the same field and line
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1),
            "Invalid coordinate 'start' on line 2 of BED file");
}

TEST(BedCtor, TruncatedGzipErrors) {
    GlobalsGuard guard;
    TempDir dir;
    std::vector<std::string> lines;
    for (int i = 0; i < 100; i++) {
        lines.push_back("chr1\t" + std::to_string(i*10) + "\t" + std::to_string(i*10 + 5));
    }
    std::string bed_fn = write_tmp_bed(dir, lines, "test.bed.gz", BEDZIP_GZIP);
    std::filesystem::resize_file(bed_fn, std::filesystem::file_size(bed_fn) / 2);

    // a half-read compressed file must not be mistaken for a short one that parsed cleanly
    EXPECT_EXIT(bedData bed(bed_fn), testing::ExitedWithCode(1), "Failed to read line");
}

TEST(BedCtor, AcceptsCarriageReturns) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t2\t8\r"});

    // htslib strips a CRLF terminator, so the stop coordinate is not read as "8\r"
    bedData bed(bed_fn);

    ASSERT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(2, bed.regions["chr1"].starts[0]);
    EXPECT_EQ(8, bed.regions["chr1"].stops[0]);
}

TEST(BedCtor, AcceptsEmptyFile) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {});

    // an empty region set is vacuously valid, and stayed so when the reader changed
    bedData bed(bed_fn);

    EXPECT_EQ(size_t(0), bed.contigs.size());
    EXPECT_EQ(0L, bed.size);
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

/* bedData::normalize *****************************************************************************/

TEST(BedNormalize, AlreadyMergedUnchanged) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}, {30, 40}});

    bed.normalize();

    EXPECT_EQ(std::vector<int>({10, 30}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({20, 40}), bed.regions["chr1"].stops);
    EXPECT_EQ(2, bed.regions["chr1"].n);
    EXPECT_EQ(20L, bed.size);
}

TEST(BedNormalize, SortsUnsorted) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{30, 40}, {10, 20}});

    // disjoint regions given out of order are ordered, not combined
    testing::internal::CaptureStderr();
    bed.normalize();
    EXPECT_NE(std::string::npos, testing::internal::GetCapturedStderr().find("were unsorted"));

    EXPECT_EQ(std::vector<int>({10, 30}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({20, 40}), bed.regions["chr1"].stops);
    EXPECT_EQ(2, bed.regions["chr1"].n);
    EXPECT_EQ(20L, bed.size);
}

TEST(BedNormalize, CombinesOverlapping) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}, {15, 30}});
    EXPECT_EQ(25L, bed.size); // add() counted [15, 20) twice

    bed.normalize();

    EXPECT_EQ(std::vector<int>({10}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({30}), bed.regions["chr1"].stops);
    EXPECT_EQ(1, bed.regions["chr1"].n);

    // the doubly-covered bases are counted once by the recomputed size
    EXPECT_EQ(20L, bed.size);
}

TEST(BedNormalize, CombinesAdjacent) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}, {20, 30}});

    // abutting regions cover a contiguous span, so they become one interval
    bed.normalize();

    EXPECT_EQ(std::vector<int>({10}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({30}), bed.regions["chr1"].stops);
    EXPECT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(20L, bed.size);
}

TEST(BedNormalize, AbsorbsNested) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 40}, {20, 30}});

    // the enclosing region must not be truncated to the nested one's stop
    bed.normalize();

    EXPECT_EQ(std::vector<int>({10}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({40}), bed.regions["chr1"].stops);
    EXPECT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(30L, bed.size);
}

TEST(BedNormalize, CollapsesDuplicates) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}, {10, 20}, {10, 20}});

    bed.normalize();

    EXPECT_EQ(std::vector<int>({10}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({20}), bed.regions["chr1"].stops);
    EXPECT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(10L, bed.size);
}

TEST(BedNormalize, SingleIntervalUnchanged) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}});

    bed.normalize();

    EXPECT_EQ(std::vector<int>({10}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({20}), bed.regions["chr1"].stops);
    EXPECT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(10L, bed.size);
}

TEST(BedNormalize, EmptyBedUnchanged) {
    GlobalsGuard guard;
    bedData bed;

    bed.normalize();

    EXPECT_EQ(size_t(0), bed.contigs.size());
    EXPECT_EQ(0L, bed.size);
}

TEST(BedNormalize, MergesEachContigSeparately) {
    GlobalsGuard guard;
    bedData bed = make_bed({{"chr1", {{10, 20}, {15, 30}}}, {"chr2", {{50, 60}}}});

    // a contig's regions never merge into another contig's, and every contig is kept
    bed.normalize();

    EXPECT_EQ(std::vector<std::string>({"chr1", "chr2"}), bed.contigs);
    EXPECT_EQ(1, bed.regions["chr1"].n);
    EXPECT_EQ(std::vector<int>({10}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({30}), bed.regions["chr1"].stops);
    EXPECT_EQ(1, bed.regions["chr2"].n);
    EXPECT_EQ(std::vector<int>({50}), bed.regions["chr2"].starts);
    EXPECT_EQ(std::vector<int>({60}), bed.regions["chr2"].stops);
    EXPECT_EQ(30L, bed.size);
}

// A merged region set is what contains() assumes, so the two must agree once merge() has run.
TEST(BedNormalize, MergedRegionsAreQueryable) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{30, 40}, {10, 25}, {20, 30}});

    testing::internal::CaptureStderr();
    bed.normalize(); // one [10, 40) region
    testing::internal::GetCapturedStderr(); // the input was unsorted, which warns

    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 26, 29, TYPE_SUB));
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 0, 5, TYPE_SUB));
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 35, 45, TYPE_SUB));
}

// The message tests load from a file rather than building in memory, since what they assert on is
// the file being named -- #47 normalizes many region sets, so a message about one of them has to
// say which.
TEST(BedNormalize, ReportsCoalescedCount) {
    GlobalsGuard guard;
    g.verbosity = 2;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t10\t20", "chr1\t15\t30", "chr1\t30\t40"});

    testing::internal::CaptureStderr();
    bedData bed(bed_fn, true); // three regions become one
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("Merged 2 overlapping or adjacent regions")) << out;
    EXPECT_NE(std::string::npos, out.find(bed_fn)) << out;

    // these regions ascend, so merging them is not grounds for reporting the file as unsorted
    EXPECT_EQ(std::string::npos, out.find("unsorted")) << out;
}

TEST(BedNormalize, WarnsWhenUnsorted) {
    GlobalsGuard guard;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t30\t40", "chr1\t10\t20"});

    // sorting is a repair, not a formality: the file is not what check() would have accepted
    testing::internal::CaptureStderr();
    bedData bed(bed_fn, true);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("[WARN")) << out;
    EXPECT_NE(std::string::npos, out.find("were unsorted")) << out;
    EXPECT_NE(std::string::npos, out.find(bed_fn)) << out;

    // the regions are sorted, and disjoint regions are not merged along the way
    EXPECT_EQ(std::vector<int>({10, 30}), bed.regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({20, 40}), bed.regions["chr1"].stops);
}

TEST(BedNormalize, SilentWhenAlreadyNormalized) {
    GlobalsGuard guard;
    g.verbosity = 2;
    TempDir dir;
    std::string bed_fn = write_tmp_bed(dir, {"chr1\t10\t20", "chr1\t30\t40"});

    // a sorted, disjoint region set is repaired in no respect, so neither message fires
    testing::internal::CaptureStderr();
    bedData bed(bed_fn, true);

    EXPECT_EQ("", testing::internal::GetCapturedStderr());
}

/* bedData::contains ******************************************************************************/

TEST(BedContains, FlippedErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                bedData bed = two_region_bed();
                bed.contains("chr1", 8, 2, TYPE_SUB);
            }, testing::ExitedWithCode(1), "Invalid region chr1:8-2 in BED contains");
}

TEST(BedContains, UnknownContigOffctg) {
    GlobalsGuard guard;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_OFFCTG, bed.contains("chrZ", 12, 15, TYPE_SUB));
}

TEST(BedContains, BeforeAllOutside) {
    GlobalsGuard guard;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 2, 5, TYPE_SUB));

    // the first region starts at 10, and a variant ending exactly there is still outside it
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 5, 10, TYPE_SUB));
}

TEST(BedContains, AfterAllOutside) {
    GlobalsGuard guard;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 45, 50, TYPE_SUB));

    // the last region ends at 40, and a variant starting exactly there is already past it
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 40, 45, TYPE_SUB));
}

TEST(BedContains, MiddleInside) {
    GlobalsGuard guard;
    bedData bed = two_region_bed();

    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 12, 15, TYPE_SUB));
    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 32, 35, TYPE_SUB));

    // a variant filling a region exactly is inside it
    EXPECT_EQ(BED_INSIDE, bed.contains("chr1", 10, 20, TYPE_SUB));
}

TEST(BedContains, InsAtRegionEndBorder) {
    GlobalsGuard guard;
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
    bedData bed = two_region_bed();

    // wholly within the [20, 30) gap
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 22, 25, TYPE_SUB));

    // flush against both gap edges, which are exclusive on the left and inclusive on the right
    EXPECT_EQ(BED_OUTSIDE, bed.contains("chr1", 20, 30, TYPE_SUB));
}

TEST(BedContains, PartialOverlapBorder) {
    GlobalsGuard guard;
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
    bedData bed = make_bed("chr1", {{10, 20}, {30, 40}, {50, 60}});

    // covers all of the middle region plus part of the outer two
    EXPECT_EQ(BED_BORDER, bed.contains("chr1", 15, 55, TYPE_SUB));
}

/* bedData::classify ******************************************************************************/

/**
 * @brief Locates a variant the way contains() does, so classify() can be called on its own.
 *
 * Mirrors the two binary searches rather than calling contains(), since a test that obtained its
 * indices from the function under test would pin nothing.
 * @param[in] bed BED holding the contig
 * @param[in] ctg Contig the variant lies on
 * @param[in] start 0-based inclusive variant start
 * @param[in] stop 0-based exclusive variant stop
 * @return Indices of the last region starting at or before start, and the first stopping at or
 *         after stop
 */
std::pair<int, int> locate(bedData & bed, const std::string & ctg, int start, int stop) {
    const std::vector<int> & starts = bed.regions[ctg].starts;
    const std::vector<int> & stops = bed.regions[ctg].stops;
    return {int(std::upper_bound(starts.begin(), starts.end(), start) - starts.begin()) - 1,
            int(std::lower_bound(stops.begin(), stops.end(), stop) - stops.begin())};
}

// The extraction is only correct if contains() and classify() cannot disagree, so every coordinate
// pair that reaches classify() through contains() must classify the same way on its own.
TEST(BedClassify, AgreesWithContains) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}, {30, 40}});

    for (int start = 0; start <= 50; start++) {
        for (int stop = start; stop <= 50; stop++) {
            for (edittype_t type : {TYPE_SUB, TYPE_INS}) {
                auto [start_idx, stop_idx] = locate(bed, "chr1", start, stop);
                EXPECT_EQ(bed.contains("chr1", start, stop, type),
                        bed.classify("chr1", start, stop, type, start_idx, stop_idx))
                        << start << "-" << stop << " type " << int(type);
            }
        }
    }
}

// A variant left of every region has start_idx -1, which the index tests would read as a partial
// overlap, so the before-all check must come first.
TEST(BedClassify, BeforeAllOutsideNotBorder) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}, {30, 40}});

    EXPECT_EQ(-1, locate(bed, "chr1", 2, 5).first);
    EXPECT_EQ(BED_OUTSIDE, bed.classify("chr1", 2, 5, TYPE_SUB, -1, 0));
}

// The mirror case: a variant right of every region has stop_idx off the end of the stop list.
TEST(BedClassify, AfterAllOutsideNotBorder) {
    GlobalsGuard guard;
    bedData bed = make_bed("chr1", {{10, 20}, {30, 40}});

    EXPECT_EQ(2, locate(bed, "chr1", 45, 50).second);
    EXPECT_EQ(BED_OUTSIDE, bed.classify("chr1", 45, 50, TYPE_SUB, 1, 2));
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

/* load_strata ************************************************************************************/

/**
 * @struct strata_result
 * @brief The names and region sets one load_strata() call produced.
 */
struct strata_result {
    std::vector<std::string> names;  ///< Stratum names, in manifest order
    std::vector<bedData> strata;     ///< Parsed stratum regions, parallel to names
};

/**
 * @brief Writes a manifest holding the given lines and loads it, returning what it parsed to.
 * @param[in] dir Temporary directory that owns the manifest, and that relative paths resolve to
 * @param[in] lines Manifest lines, written verbatim in order
 * @return The loaded stratum names and region sets
 */
strata_result load_manifest(const TempDir & dir, const std::vector<std::string> & lines) {
    strata_result result;
    load_strata(write_tmp_text(dir, lines, "strata.tsv"), result.names, result.strata);
    return result;
}

TEST(LoadStrata, ManifestOpenFailErrors) {
    GlobalsGuard guard;
    TempDir dir;
    std::vector<std::string> names;
    std::vector<bedData> strata;

    EXPECT_EXIT(load_strata(dir.path("absent.tsv"), names, strata), testing::ExitedWithCode(1),
            "Failed to open stratification manifest");
}

TEST(LoadStrata, ResolvesRelativeToManifestDir) {
    GlobalsGuard guard;
    TempDir dir;
    std::filesystem::create_directories(dir.path("sets"));
    write_tmp_bed(dir, {"chr1\t10\t20"}, "sets/ancestry.bed");

    // "sets/ancestry.bed" names nothing beneath the working directory, so only resolution against
    // the manifest's own directory can find it -- which is what the GIAB manifests require
    strata_result r = load_manifest(dir, {"ancestry_AFR\tsets/ancestry.bed"});

    EXPECT_EQ(std::vector<std::string>({"ancestry_AFR"}), r.names);
    ASSERT_EQ(size_t(1), r.strata.size());
    EXPECT_EQ(std::vector<std::string>({"chr1"}), r.strata[0].contigs);
    EXPECT_EQ(10L, r.strata[0].size);
}

TEST(LoadStrata, AcceptsAbsolutePaths) {
    GlobalsGuard guard;
    TempDir dir;
    const std::string bed_fn = write_tmp_bed(dir, {"chr1\t10\t20"}, "absolute.bed");
    ASSERT_TRUE(std::filesystem::path(bed_fn).is_absolute());

    // an absolute path is opened as written, rather than being appended to the manifest's directory
    strata_result r = load_manifest(dir, {"refseq_cds\t" + bed_fn});

    EXPECT_EQ(std::vector<std::string>({"refseq_cds"}), r.names);
    ASSERT_EQ(size_t(1), r.strata.size());
    EXPECT_EQ(10L, r.strata[0].size);
}

TEST(LoadStrata, SkipsCommentsAndBlankLines) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t10\t20"}, "a.bed");

    strata_result r = load_manifest(dir,
            {"# GRCh38 stratifications", "", "a\ta.bed", "", "#trailing comment"});

    EXPECT_EQ(std::vector<std::string>({"a"}), r.names);
    EXPECT_EQ(size_t(1), r.strata.size());
}

TEST(LoadStrata, IgnoresExtraColumns) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t10\t20"}, "a.bed");

    strata_result r = load_manifest(dir, {"a\ta.bed\tsome\tother\tcolumns"});

    EXPECT_EQ(std::vector<std::string>({"a"}), r.names);
    ASSERT_EQ(size_t(1), r.strata.size());
    EXPECT_EQ(10L, r.strata[0].size);
}

TEST(LoadStrata, PreservesManifestOrder) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t10\t20"}, "z.bed");
    write_tmp_bed(dir, {"chr1\t30\t50"}, "a.bed");

    // load order is manifest order, since it is also the output row order the strata will take
    strata_result r = load_manifest(dir, {"zzz\tz.bed", "aaa\ta.bed"});

    EXPECT_EQ(std::vector<std::string>({"zzz", "aaa"}), r.names);
    ASSERT_EQ(size_t(2), r.strata.size());
    EXPECT_EQ(10L, r.strata[0].size);
    EXPECT_EQ(20L, r.strata[1].size);
}

TEST(LoadStrata, NormalizesRegionsAtLoad) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t30\t40", "chr1\t10\t35"}, "messy.bed");

    // a third-party region set is normalized rather than rejected, so an unsorted and overlapping
    // file loads as one merged interval instead of failing the check the same file would fail as
    // an evaluation BED; the sort is still reported, since the file is not what check() accepts
    testing::internal::CaptureStderr();
    strata_result r = load_manifest(dir, {"messy\tmessy.bed"});
    std::string out = testing::internal::GetCapturedStderr();

    ASSERT_EQ(size_t(1), r.strata.size());
    EXPECT_EQ(std::vector<int>({10}), r.strata[0].regions["chr1"].starts);
    EXPECT_EQ(std::vector<int>({40}), r.strata[0].regions["chr1"].stops);
    EXPECT_EQ(30L, r.strata[0].size);
    EXPECT_NE(std::string::npos, out.find("were unsorted, and have been sorted")) << out;
}

// The strata are normalized, not exempted from validation: a malformation sorting and merging
// cannot repair is as fatal in a stratum BED as in an evaluation BED.
TEST(LoadStrata, StratumBedStillChecked) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t8\t2"}, "flipped.bed");

    EXPECT_EXIT(load_manifest(dir, {"flipped\tflipped.bed"}), testing::ExitedWithCode(1),
            "BED region chr1:8-2 stop precedes start");
}

TEST(LoadStrata, MissingFieldErrors) {
    GlobalsGuard guard;
    TempDir dir;

    // line numbers count skipped lines too, so the reported line is the one in the file
    EXPECT_EXIT(load_manifest(dir, {"# comment", "name_with_no_path"}),
            testing::ExitedWithCode(1),
            "Line 2 of stratification manifest .* does not name both a stratum and a BED file");
}

TEST(LoadStrata, EmptyPathFieldErrors) {
    GlobalsGuard guard;
    TempDir dir;

    // a trailing tab supplies a second field that names nothing, which is not a BED path
    EXPECT_EXIT(load_manifest(dir, {"name\t"}), testing::ExitedWithCode(1),
            "Line 1 of stratification manifest .* does not name both a stratum and a BED file");
}

TEST(LoadStrata, DuplicateNameErrors) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t10\t20"}, "a.bed");

    // two rows sharing a name would make every stratified output row ambiguous
    EXPECT_EXIT(load_manifest(dir, {"dup\ta.bed", "dup\ta.bed"}), testing::ExitedWithCode(1),
            "Duplicate stratum name 'dup' on line 2 of stratification manifest");
}

TEST(LoadStrata, ReservedStarNameErrors) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t10\t20"}, "a.bed");

    EXPECT_EXIT(load_manifest(dir, {"*\ta.bed"}), testing::ExitedWithCode(1),
            "Stratum name '\\*' on line 1 of stratification manifest .* is reserved");
}

TEST(LoadStrata, StratumBedOpenFailErrors) {
    GlobalsGuard guard;
    TempDir dir;

    // both paths are named: the manifest-relative one identifies the offending row, and the
    // resolved one shows where it was looked for, which a relative path alone cannot
    EXPECT_EXIT(load_manifest(dir, {"absent\tsets/absent.bed"}), testing::ExitedWithCode(1),
            "Failed to open BED file 'sets/absent.bed' for stratum 'absent' on line 1 of"
            " stratification manifest .* resolved to '/.*/sets/absent.bed'");
}

TEST(LoadStrata, NoStrataWarns) {
    GlobalsGuard guard;
    TempDir dir;

    // a manifest naming nothing is not an error; the run proceeds as though -st were absent
    testing::internal::CaptureStderr();
    strata_result r = load_manifest(dir, {"# every line is a comment", ""});
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_EQ(size_t(0), r.names.size());
    EXPECT_EQ(size_t(0), r.strata.size());
    EXPECT_NE(std::string::npos, out.find("[WARN")) << out;
    EXPECT_NE(std::string::npos, out.find("names no region sets")) << out;
}

TEST(LoadStrata, ReplacesPriorContents) {
    GlobalsGuard guard;
    TempDir dir;
    write_tmp_bed(dir, {"chr1\t10\t20"}, "a.bed");
    write_tmp_bed(dir, {"chr2\t10\t20"}, "b.bed");
    std::vector<std::string> names;
    std::vector<bedData> strata;

    load_strata(write_tmp_text(dir, {"a\ta.bed"}, "first.tsv"), names, strata);
    load_strata(write_tmp_text(dir, {"b\tb.bed"}, "second.tsv"), names, strata);

    // a second load replaces the strata rather than appending to them
    EXPECT_EQ(std::vector<std::string>({"b"}), names);
    ASSERT_EQ(size_t(1), strata.size());
    EXPECT_EQ(std::vector<std::string>({"chr2"}), strata[0].contigs);
}

/* check_strata_contigs ***************************************************************************/

TEST(CheckStrataContigs, NoStrataSilent) {
    GlobalsGuard guard;

    testing::internal::CaptureStderr();
    check_strata_contigs(two_contig_ref());

    EXPECT_EQ("", testing::internal::GetCapturedStderr());
}

TEST(CheckStrataContigs, SharedContigSilent) {
    GlobalsGuard guard;
    g.strat_names = {"shared"};
    g.strata = {make_bed({{"chr2", {{0, 10}}}, {"chrZ", {{0, 10}}}})};

    // one shared contig is enough; the region set may name contigs the reference does not have
    testing::internal::CaptureStderr();
    check_strata_contigs(two_contig_ref());

    EXPECT_EQ("", testing::internal::GetCapturedStderr());
}

TEST(CheckStrataContigs, ZeroOverlapWarnsByName) {
    GlobalsGuard guard;
    g.ref_fasta_fn = "ref.fasta";
    g.strat_names = {"matching", "offctg"};
    g.strata = {make_bed("chr1", {{0, 10}}), make_bed("1", {{0, 10}})};

    testing::internal::CaptureStderr();
    check_strata_contigs(two_contig_ref());
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("[WARN")) << out;
    EXPECT_NE(std::string::npos, out.find("Stratification region set 'offctg' shares no contig"
            " with reference FASTA 'ref.fasta'")) << out;
    EXPECT_EQ(std::string::npos, out.find("'matching'")) << out;
}

TEST(CheckStrataContigs, AllZeroOverlapNamesLikelyCause) {
    GlobalsGuard guard;
    g.ref_fasta_fn = "ref.fasta";
    g.strat_tsv_fn = "strata.tsv";
    g.strat_names = {"one", "two"};
    g.strata = {make_bed("1", {{0, 10}}), make_bed("2", {{0, 10}})};

    // every set failing at once is an assembly or contig-naming mismatch rather than 181
    // independent accidents, so the cause is named instead of each set
    testing::internal::CaptureStderr();
    check_strata_contigs(two_contig_ref());
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("No stratification region set in 'strata.tsv' shares a"
            " contig with reference FASTA 'ref.fasta'")) << out;
    EXPECT_NE(std::string::npos, out.find("contig naming")) << out;
    EXPECT_EQ(std::string::npos, out.find("Stratification region set 'one'")) << out;
}

/* intersect_contigs ******************************************************************************/

// A contig outside the BED survives so that its variants remain reportable. It is already empty:
// parse_variants() discards every variant on it as BED_OFFCTG, so keeping it costs an entry.
TEST(IntersectContigs, BedKeepsContigOutsideBed) {
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

    // chr2 keeps its entry in contigs and in every field parallel to it
    EXPECT_EQ(std::vector<std::string>({"chr1", "chr2"}), query->contigs);
    EXPECT_EQ(std::vector<int>({12, 12}), query->lengths);
    EXPECT_EQ(std::vector< std::set<int> >({{2}, {2}}), query->observed_ploidies);
    EXPECT_EQ(size_t(1), query->variants[HAP1].count("chr2"));
    EXPECT_EQ(size_t(1), query->variants[HAP2].count("chr2"));

    // set_var_record() anchors an INS/DEL in the reference sequence, so it must survive too
    EXPECT_EQ(size_t(2), ref->fasta.size());
    EXPECT_EQ("TTTTTTTTTTTT", ref->fasta.at("chr2"));
}

// The BED contig list, not the truth VCF, decides which reference sequences are needed, and neither
// decides what may be freed: a contig absent from both is still reachable through a retained call.
TEST(IntersectContigs, NobedKeepsRefContigAbsentFromTruth) {
    GlobalsGuard guard;
    g.bed_exists = false;
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
    std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    testing::internal::GetCapturedStderr();

    EXPECT_EQ(size_t(2), ref->fasta.size());
    EXPECT_EQ("TTTTTTTTTTTT", ref->fasta.at("chr2"));
}

// Retaining a contig means its reference sequence is now required, so a missing one is an error
// rather than a silent drop. Checked for both callsets, and whether or not a BED was supplied.
TEST(IntersectContigs, QueryContigMissingFromFastaErrors) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{0, 10}});
    EXPECT_EXIT({
                std::shared_ptr<variantData> query =
                        make_variantData(QUERY, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
                std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
                std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
                intersect_contigs(query, truth, ref);
            }, testing::ExitedWithCode(1),
            "Contig 'chr2' found in QUERY VCF but not reference FASTA");
}

TEST(IntersectContigs, BedTruthContigMissingFromFastaErrors) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{0, 10}});
    EXPECT_EXIT({
                std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
                std::shared_ptr<variantData> truth =
                        make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
                std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
                intersect_contigs(query, truth, ref);
            }, testing::ExitedWithCode(1),
            "Contig 'chr2' found in TRUTH VCF but not reference FASTA");
}

// Nothing is evaluated on a contig outside the BED, so there are no false positives or false
// negatives there to warn about.
TEST(IntersectContigs, ContigOutsideBedDoesNotWarnOnMissingCounterpart) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{0, 10}});
    std::shared_ptr<variantData> query =
            make_variantData(QUERY, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
    std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_EQ(std::string::npos, out.find("Contig 'chr2'")) << out;
}

// Ploidies are recorded before parse_variants() applies the BED filter, so a contig outside the BED
// carries observed ploidies that no evaluation ever consults. Disagreement there is not actionable.
TEST(IntersectContigs, ContigOutsideBedDoesNotWarnOnPloidy) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{0, 10}});
    std::shared_ptr<variantData> query =
            make_variantData(QUERY, {"chr1", "chr2"}, {12, 12}, {{2}, {1}});
    std::shared_ptr<variantData> truth =
            make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_EQ(std::string::npos, out.find("contig 'chr2' has ploid")) << out;
}

// A truth contig outside the BED has no query counterpart to look up, because the injection block
// below iterates the BED contigs only. Indexing observed_ploidies with the result of a failed
// std::find would read past the end of the vector.
TEST(IntersectContigs, TruthOnlyContigOutsideBedIsSafe) {
    GlobalsGuard guard;
    g.bed_exists = true;
    g.bed = make_bed("chr1", {{0, 10}});
    std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
    std::shared_ptr<variantData> truth =
            make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {1}});
    std::shared_ptr<fastaData> ref = two_contig_ref();

    testing::internal::CaptureStderr();
    intersect_contigs(query, truth, ref);
    std::string out = testing::internal::GetCapturedStderr();

    // chr2 stays truth-only: the BED does not name it, so nothing injects a query counterpart
    EXPECT_EQ(std::vector<std::string>({"chr1"}), query->contigs);
    EXPECT_EQ(std::string::npos, out.find("contig 'chr2' has ploid")) << out;
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

TEST(IntersectContigs, NobedTruthContigMissingFromFastaErrors) {
    GlobalsGuard guard;
    g.bed_exists = false;
    EXPECT_EXIT({
                std::shared_ptr<variantData> query = make_variantData(QUERY, {"chr1"}, {12}, {{2}});
                std::shared_ptr<variantData> truth =
                        make_variantData(TRUTH, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
                std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
                intersect_contigs(query, truth, ref);
            }, testing::ExitedWithCode(1),
            "Contig 'chr2' found in TRUTH VCF but not reference FASTA");
}

// Without a BED the missing sequence was previously reached as an unhandled std::out_of_range from
// lengths.at() while injecting the contig into the truth, rather than as a reported error.
TEST(IntersectContigs, NobedQueryContigMissingFromFastaErrors) {
    GlobalsGuard guard;
    g.bed_exists = false;
    EXPECT_EXIT({
                std::shared_ptr<variantData> query =
                        make_variantData(QUERY, {"chr1", "chr2"}, {12, 12}, {{2}, {2}});
                std::shared_ptr<variantData> truth = make_variantData(TRUTH, {"chr1"}, {12}, {{2}});
                std::shared_ptr<fastaData> ref = make_fasta("chr1", "ACGTACGTACGT");
                intersect_contigs(query, truth, ref);
            }, testing::ExitedWithCode(1),
            "Contig 'chr2' found in QUERY VCF but not reference FASTA");
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
