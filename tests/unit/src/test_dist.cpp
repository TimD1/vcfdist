#include <string>
#include <unordered_set>

#include "gtest/gtest.h"
#include "../../../src/dist.h"
#include "../../../src/globals.h"

namespace {

// contains is instantiated below over element types dist.cpp never uses, so these tests only link
// while its definition lives in dist.h rather than dist.cpp.

TEST(Contains, SetInt) {
    std::unordered_set<int> hash_set = {-7, 0, 42};
    EXPECT_TRUE(contains(hash_set, -7));
    EXPECT_TRUE(contains(hash_set, 0));
    EXPECT_TRUE(contains(hash_set, 42));
    EXPECT_FALSE(contains(hash_set, 7));
    EXPECT_FALSE(contains(hash_set, -42));
}

TEST(Contains, SetIntEmpty) {
    std::unordered_set<int> hash_set;
    EXPECT_FALSE(contains(hash_set, 0));
}

TEST(Contains, SetString) {
    std::unordered_set<std::string> hash_set = {"", "ACGT"};
    EXPECT_TRUE(contains(hash_set, std::string("")));
    EXPECT_TRUE(contains(hash_set, std::string("ACGT")));
    // matches are exact: lowercase, prefixes, and extensions are all absent
    EXPECT_FALSE(contains(hash_set, std::string("acgt")));
    EXPECT_FALSE(contains(hash_set, std::string("ACG")));
    EXPECT_FALSE(contains(hash_set, std::string("ACGTA")));
}

TEST(TestNG50, TestNG50Calc) {
    std::vector<int> phase_blocks1 = {25, 50};
    EXPECT_EQ(50, calc_ng50(phase_blocks1, 100));

    std::vector<int> phase_blocks2 = {1, 25, 25, 50};
    EXPECT_EQ(25, calc_ng50(phase_blocks2, 101));

    std::vector<int> phase_blocks3 = {90};
    EXPECT_EQ(90, calc_ng50(phase_blocks3, 100));

    std::vector<int> phase_blocks4 = {1, 10, 20};
    EXPECT_EQ(0, calc_ng50(phase_blocks4, 100));
}

} // namespace
