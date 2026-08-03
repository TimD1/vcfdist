/**
 * @file test_globals.cpp
 * @brief Unit tests for globals.cpp: parent_path, create_directory, init_timers, string tables.
 * @note `parse_args` is deliberately not covered here; its cases are gated on pending changes to
 *       the command-line flag surface.
 */
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/defs.h"
#include "../../../src/globals.h"
#include "../../../src/timer.h"
#include "test_helpers.h"

namespace {

/** @brief Returns the names of the given timers, in order. */
std::vector<std::string> timer_names(std::vector<timer> & timers) {
    std::vector<std::string> names;
    for (timer & t : timers) names.push_back(t.get_name());
    return names;
}

/* parent_path ************************************************************************************/

TEST(ParentPath, Nested) {
    EXPECT_EQ("a/b/", parent_path("a/b/c"));
}

TEST(ParentPath, SingleDir) {
    EXPECT_EQ("out/", parent_path("out/prefix"));
}

TEST(ParentPath, BareFilename) {
    // no separator anywhere, so there is no parent to return
    EXPECT_EQ("", parent_path("prefix"));
}

TEST(ParentPath, TrailingSlash) {
    // the trailing separator is itself the last one found, so the path is returned unchanged
    EXPECT_EQ("a/b/", parent_path("a/b/"));
}

TEST(ParentPath, Root) {
    EXPECT_EQ("/", parent_path("/"));
}

TEST(ParentPath, AbsoluteFile) {
    EXPECT_EQ("/", parent_path("/prefix"));
}

TEST(ParentPath, Empty) {
    EXPECT_EQ("", parent_path(""));
}

TEST(ParentPath, DotSlash) {
    EXPECT_EQ("./", parent_path("./prefix"));
}

TEST(ParentPath, DotDot) {
    EXPECT_EQ("../out/", parent_path("../out/prefix"));
}

/* create_directory *******************************************************************************/

TEST(CreateDirectory, SingleDir) {
    TempDir dir;
    create_directory(dir.path("a/b"));

    // only components followed by a separator are created, so the final "b" is not
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("a")));
    EXPECT_FALSE(std::filesystem::exists(dir.path("a/b")));
}

TEST(CreateDirectory, NestedTrailingSlash) {
    TempDir dir;
    create_directory(dir.path("x/y/"));

    // the trailing separator makes "y" an intermediate component, so every level is created
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("x")));
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("x/y")));
}

TEST(CreateDirectory, AlreadyExistsOk) {
    TempDir dir;

    // every component of the first call already exists on the second, and EEXIST is tolerated
    create_directory(dir.path("dup/"));
    create_directory(dir.path("dup/"));
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("dup")));
}

TEST(CreateDirectory, AbsoluteSkipsRoot) {
    // the scan starts one character past the front of the path, so a leading separator is never
    // treated as one: a single top-level component yields no mkdir call at all
    create_directory("/vcfdist_create_directory_should_be_a_noop");
    EXPECT_FALSE(std::filesystem::exists("/vcfdist_create_directory_should_be_a_noop"));

    // and a deeper absolute path succeeds without ever attempting to create the root
    TempDir dir;
    create_directory(dir.path("abs/leaf"));
    EXPECT_TRUE(std::filesystem::is_directory(dir.path("abs")));
}

TEST(CreateDirectory, MkdirFailureErrors) {
    GlobalsGuard guard;
    TempDir dir;

    // a regular file where a directory is needed makes mkdir fail with ENOTDIR, not EEXIST
    std::ofstream(dir.path("file")) << "not a directory\n";
    ASSERT_TRUE(std::filesystem::is_regular_file(dir.path("file")));

    EXPECT_EXIT(create_directory(dir.path("file/sub/")), testing::ExitedWithCode(1),
            "Unable to create directory");
}

TEST(CreateDirectory, NoSlashNoop) {
    const std::string relative = "vcfdist_create_directory_no_slash_noop";
    ASSERT_FALSE(std::filesystem::exists(relative));

    // with no separator there is nothing to iterate over, so no directory is created
    create_directory(relative);
    EXPECT_FALSE(std::filesystem::exists(relative));
}

/* init_timers ************************************************************************************/

TEST(InitTimers, Populates) {
    GlobalsGuard guard;
    g.timers.clear();

    g.init_timers(timer_strs);

    ASSERT_EQ(size_t(TIME_TOTAL+1), g.timers.size());
    EXPECT_EQ(timer_strs, timer_names(g.timers));
    EXPECT_EQ("reading", g.timers[TIME_READ].get_name());
    EXPECT_EQ("total", g.timers[TIME_TOTAL].get_name());
}

TEST(InitTimers, EmptyInput) {
    GlobalsGuard guard;
    g.timers.clear();

    g.init_timers({});

    EXPECT_TRUE(g.timers.empty());
}

TEST(InitTimers, AppendsNotClears) {
    GlobalsGuard guard;
    g.timers.clear();

    // each call pushes onto the existing vector rather than replacing it
    g.init_timers({"first"});
    g.init_timers({"second", "third"});

    ASSERT_EQ(size_t(3), g.timers.size());
    EXPECT_EQ(std::vector<std::string>({"first", "second", "third"}), timer_names(g.timers));
}

TEST(InitTimers, WritesThisNotGlobal) {
    GlobalsGuard guard;
    g.timers.clear();

    // the timers land on the instance the method was called on, not on the global `g`
    Globals local;
    local.init_timers({"local"});

    ASSERT_EQ(size_t(1), local.timers.size());
    EXPECT_EQ("local", local.timers[0].get_name());
    EXPECT_TRUE(g.timers.empty());
}

/* String lookup tables ***************************************************************************/

TEST(StringTables, SizesMatchCount) {
    // these five tables are sized exactly to their governing count constant
    EXPECT_EQ(size_t(TYPES), type_strs.size());
    EXPECT_EQ(size_t(VARTYPES), vartype_strs.size());
    EXPECT_EQ(size_t(ERRTYPES), error_strs.size());
    EXPECT_EQ(size_t(SWITCHTYPES), switch_strs.size());
    EXPECT_EQ(size_t(CALLSETS), callset_strs.size());
}

TEST(StringTables, SizesWithSentinel) {
    // these four are one longer than their count constant, because a sentinel value is also a
    // valid subscript; shortening any of them to its count would read out of bounds
    EXPECT_EQ(size_t(AC_ERRTYPES+1), ac_strs.size());
    EXPECT_EQ(size_t(PHASES+1), phase_strs.size());
    EXPECT_EQ(".", ac_strs[AC_UNKNOWN]);
    EXPECT_EQ(".", phase_strs[PHASE_NONE]);

    // gt_strs and timer_strs have no count constant, so the highest valid index bounds them
    EXPECT_EQ(size_t(GT_OTHER+1), gt_strs.size());
    EXPECT_EQ(size_t(TIME_TOTAL+1), timer_strs.size());
    EXPECT_EQ("M|N", gt_strs[GT_OTHER]);
    EXPECT_EQ("total", timer_strs[TIME_TOTAL]);
}

TEST(StringTables, IndexMapping) {
    EXPECT_EQ("TP", error_strs[ERRTYPE_TP]);
    EXPECT_EQ("1|1", gt_strs[GT_ALT1_ALT1]);
    EXPECT_EQ("SNP", type_strs[TYPE_SUB]);
    EXPECT_EQ("QUERY", callset_strs[QUERY]);
    EXPECT_EQ("TRUTH", callset_strs[TRUTH]);
    EXPECT_EQ("SV", vartype_strs[VARTYPE_SV]);
    EXPECT_EQ("SWITCH", switch_strs[SWITCHTYPE_SWITCH]);
    EXPECT_EQ("INSIDE ", region_strs[BED_INSIDE]);
}

TEST(StringTables, AliasedIndices) {
    // several constant pairs deliberately share a subscript
    EXPECT_EQ(TYPE_REF, TYPE_ALL);
    EXPECT_EQ(TYPE_CPX, TYPE_INDEL);
    EXPECT_EQ(ERRTYPE_UN, ERRTYPE_NE);
    EXPECT_EQ(REF, TRUTH);

    // type_strs and type_strs2 are parallel but not interchangeable: they disagree at the
    // aliased subscripts, since type_strs2 names the aggregation class instead of the variant type
    ASSERT_EQ(type_strs.size(), type_strs2.size());
    EXPECT_EQ("REF", type_strs[TYPE_REF]);
    EXPECT_EQ("ALL", type_strs2[TYPE_ALL]);
    EXPECT_EQ("CPX", type_strs[TYPE_CPX]);
    EXPECT_EQ("INDEL", type_strs2[TYPE_INDEL]);

    // the aliases agree elsewhere, so only indices 0 and 4 differ
    EXPECT_EQ(type_strs[TYPE_SUB], type_strs2[TYPE_SUB]);
    EXPECT_EQ(type_strs[TYPE_INS], type_strs2[TYPE_INS]);
    EXPECT_EQ(type_strs[TYPE_DEL], type_strs2[TYPE_DEL]);

    // a single string serves both unknown and not-evaluated error types
    EXPECT_EQ("??", error_strs[ERRTYPE_UN]);
}

TEST(StringTables, AcSparse) {
    // only the two transitions that change the allele count while still counting as a TP carry a
    // non-"." string; every other allele-count error type prints as "."
    EXPECT_EQ("+", ac_strs[AC_ERR_1_TO_2]);
    EXPECT_EQ("-", ac_strs[AC_ERR_2_TO_1]);
    for (size_t i = 0; i < ac_strs.size(); i++) {
        if (i == AC_ERR_1_TO_2 || i == AC_ERR_2_TO_1) continue;
        EXPECT_EQ(".", ac_strs[i]) << "ac_strs[" << i << "]";
    }
}

TEST(StringTables, RegionPadded) {
    // entries are space-padded to a common width, which downstream column output relies on; the
    // padding surfaces as a trailing space in the LOCATION column of query.tsv and truth.tsv
    EXPECT_EQ("OUTSIDE", region_strs[BED_OUTSIDE]);
    EXPECT_EQ("INSIDE ", region_strs[BED_INSIDE]);
    EXPECT_EQ("BORDER ", region_strs[BED_BORDER]);
    EXPECT_EQ("OFF CTG", region_strs[BED_OFFCTG]);

    ASSERT_FALSE(region_strs.empty());
    for (const std::string & s : region_strs) {
        EXPECT_EQ(region_strs[0].size(), s.size()) << "unpadded entry '" << s << "'";
    }
}

} // namespace
