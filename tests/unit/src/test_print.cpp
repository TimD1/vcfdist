/**
 * @file test_print.cpp
 * @brief Unit tests for print.cpp: qscore, get_ptr_repr, color wrappers, write_params.
 * @note compute_pr_f1 is not covered here; its extraction is tracked separately in issue #94.
 */
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/dist.h"
#include "../../../src/globals.h"
#include "../../../src/print.h"
#include "test_helpers.h"

namespace {

/**
 * @brief Reads an entire file into a string.
 * @param[in] fn Path of the file to read
 * @return File contents, or an empty string if the file cannot be opened
 */
std::string read_file(const std::string & fn) {
    std::ifstream in(fn);
    if (!in.is_open()) return "";
    std::ostringstream ostr;
    ostr << in.rdbuf();
    return ostr.str();
}

/**
 * @brief Builds a single-entry pointer map for get_ptr_repr.
 * @param[in] cell Cell to look up
 * @param[in] prev Predecessor cell the map returns
 * @return Map containing exactly the one cell-to-predecessor entry
 */
std::unordered_map<idx4,idx4> one_ptr(const idx4 & cell, const idx4 & prev) {
    std::unordered_map<idx4,idx4> ptrs;
    ptrs[cell] = prev;
    return ptrs;
}

/* qscore *****************************************************************************************/

TEST(Qscore, Prob1) {
    EXPECT_FLOAT_EQ(0.0f, qscore(1.0));
}

TEST(Qscore, Prob0p1) {
    EXPECT_FLOAT_EQ(10.0f, qscore(0.1));
}

TEST(Qscore, Prob0p01) {
    EXPECT_FLOAT_EQ(20.0f, qscore(0.01));
}

TEST(Qscore, Prob0p001) {
    EXPECT_FLOAT_EQ(30.0f, qscore(0.001));
}

TEST(Qscore, ClampHigh) {
    // -10*log10(1e-11) is 110, above the upper clamp; the cap is 100, not the Phred-60 convention
    EXPECT_FLOAT_EQ(100.0f, qscore(1e-11));
}

TEST(Qscore, Prob0) {
    // log10(0) is -inf, so the score is +inf before the upper clamp reduces it to 100
    EXPECT_FLOAT_EQ(100.0f, qscore(0.0));
}

TEST(Qscore, NegativeInput) {
    // log10 of a negative is NaN; std::max(0.0, NaN) returns its first argument, so the result is 0
    float q = qscore(-1.0);
    EXPECT_FALSE(std::isnan(q));
    EXPECT_FLOAT_EQ(0.0f, q);
}

TEST(Qscore, ProbGt1) {
    // -10*log10(2) is negative, so the lower clamp returns 0
    EXPECT_FLOAT_EQ(0.0f, qscore(2.0));
}

TEST(Qscore, Rounding) {
    // the return type is float, so fractional scores are preserved rather than rounded to integers
    float q = qscore(0.5);
    EXPECT_NEAR(3.0103f, q, 1e-4);
    EXPECT_NE(0.0f, q - std::floor(q));
}

/* get_ptr_repr ***********************************************************************************/

TEST(GetPtrRepr, NotFound) {
    idx4 cell(0, 0, 1, 1);
    std::unordered_map<idx4,idx4> ptrs;
    EXPECT_EQ("  .", get_ptr_repr(cell, ptrs));
}

TEST(GetPtrRepr, Up) {
    idx4 cell(0, 0, 3, 5);
    idx4 prev(0, 0, 2, 5);
    EXPECT_EQ("  |", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, Left) {
    idx4 cell(0, 0, 3, 5);
    idx4 prev(0, 0, 3, 4);
    EXPECT_EQ("  _", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, Diag) {
    idx4 cell(0, 0, 3, 5);
    idx4 prev(0, 0, 2, 4);
    EXPECT_EQ("  \\", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, InvalidSameMatrix) {
    // same node pair, but the predecessor is neither adjacent nor diagonal
    idx4 cell(0, 0, 5, 5);
    idx4 prev(0, 0, 1, 1);
    EXPECT_EQ(" ?1", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, QueryNode) {
    // entering query node 8 from query node 7 prints the predecessor id, zero-padded to 2 digits
    idx4 cell(8, 0, 0, 4);
    idx4 prev(7, 0, 2, 4);
    EXPECT_EQ("^07", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, TruthNode) {
    // entering truth node 8 from truth node 3 prints the predecessor id, zero-padded to 2 digits
    idx4 cell(0, 8, 3, 0);
    idx4 prev(0, 3, 3, 2);
    EXPECT_EQ("<03", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, InvalidOther) {
    // the query node changed, but not at the start of the node, so no node-id branch applies
    idx4 cell(1, 0, 5, 0);
    idx4 prev(0, 0, 4, 0);
    EXPECT_EQ(" ?2", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, WidthConsistency) {
    // every branch must return exactly 3 characters for node ids below 100 so the grid stays
    // aligned, including the two invalid markers (issue #73)
    std::vector< std::pair<idx4, idx4> > cases = {
        {idx4(0, 0, 3, 5), idx4(0, 0, 2, 5)}, // up
        {idx4(0, 0, 3, 5), idx4(0, 0, 3, 4)}, // left
        {idx4(0, 0, 3, 5), idx4(0, 0, 2, 4)}, // diagonal
        {idx4(0, 0, 5, 5), idx4(0, 0, 1, 1)}, // invalid, same matrix
        {idx4(8, 0, 0, 4), idx4(7, 0, 2, 4)}, // new query node, two digits
        {idx4(1, 0, 0, 4), idx4(0, 0, 2, 4)}, // new query node, one digit
        {idx4(0, 8, 3, 0), idx4(0, 3, 3, 2)}, // new truth node, two digits
        {idx4(0, 1, 3, 0), idx4(0, 0, 3, 2)}, // new truth node, one digit
        {idx4(1, 0, 5, 0), idx4(0, 0, 4, 0)}, // invalid, different matrix
    };
    for (size_t i = 0; i < cases.size(); i++) {
        std::string repr = get_ptr_repr(cases[i].first, one_ptr(cases[i].first, cases[i].second));
        EXPECT_EQ(size_t(3), repr.size()) << "case " << i << ": '" << repr << "'";
    }

    // the not-found branch takes no pointer map entry
    std::unordered_map<idx4,idx4> empty;
    EXPECT_EQ(size_t(3), get_ptr_repr(idx4(0, 0, 1, 1), empty).size());
}

TEST(GetPtrRepr, QueryNodeWide) {
    // std::setw(2) pads but never truncates, so a node id of 100 or more widens the cell to 4
    // characters and misaligns the debug grid; this pins current output rather than correctness
    idx4 cell(124, 0, 0, 4);
    idx4 prev(123, 0, 2, 4);
    std::string repr = get_ptr_repr(cell, one_ptr(cell, prev));
    EXPECT_EQ("^123", repr);
    EXPECT_EQ(size_t(4), repr.size());
}

/* color wrappers *********************************************************************************/

TEST(Color, Green) {
    // the wrappers emit escape codes unconditionally; the isatty check lives in the COLOR_* macros
    EXPECT_EQ("\033[32m7\033[0m", GREEN(7));
}

TEST(Color, Red) {
    EXPECT_EQ("\033[31mA\033[0m", RED('A'));
}

TEST(Color, Blue) {
    EXPECT_EQ("\033[34mACGT\033[0m", BLUE(std::string("ACGT")));
}

TEST(Color, Yellow) {
    EXPECT_EQ("\033[33m-3\033[0m", YELLOW(-3));
}

TEST(Color, Purple) {
    EXPECT_EQ("\033[35mchr1\033[0m", PURPLE(std::string("chr1")));
}

TEST(Color, GreenStrEmpty) {
    // an empty payload still yields the prefix and reset suffix
    EXPECT_EQ("\033[32m\033[0m", GREEN(std::string("")));
}

/* write_params ***********************************************************************************/

TEST(WriteParams, FiltersJoin) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";
    g.filters = {"A", "B", "C"};

    write_params();

    std::string contents = read_file(g.out_prefix + "parameters.tsv");
    EXPECT_NE(std::string::npos, contents.find("\nfilters\tA,B,C\n")) << contents;
}

TEST(WriteParams, FiltersEmpty) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";
    g.filters.clear();

    // guards the out-of-bounds g.filters[0] read that an unguarded join would perform
    write_params();

    std::string contents = read_file(g.out_prefix + "parameters.tsv");
    EXPECT_NE(std::string::npos, contents.find("\nfilters\t\n")) << contents;
}

TEST(WriteParams, FiltersSingle) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";
    g.filters = {"PASS"};

    write_params();

    std::string contents = read_file(g.out_prefix + "parameters.tsv");
    EXPECT_NE(std::string::npos, contents.find("\nfilters\tPASS\n")) << contents;
    EXPECT_EQ(std::string::npos, contents.find("PASS,"));
}

TEST(WriteParams, FopenFail) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path("missing/");

    EXPECT_EXIT(write_params(), testing::ExitedWithCode(1), "Failed to open parameters TSV file");
}

} // namespace
