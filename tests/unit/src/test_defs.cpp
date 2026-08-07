/**
 * @file test_defs.cpp
 * @brief Unit tests for the enum-keyed containers in defs.h.
 */
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "../../../src/defs.h"

namespace {

// a throwaway family, so these tests do not move when a real family gains an enumerator
enum class fruit_t : int8_t { APPLE = 0, PEAR = 1, PLUM = 2 };
constexpr std::size_t FRUIT_SLOTS = 3;

enum class veg_t : int8_t { KALE = 0 };

// detects whether an array of type A accepts a subscript of type K, without instantiating it
template <typename A, typename K, typename = void>
struct has_subscript : std::false_type {};
template <typename A, typename K>
struct has_subscript<A, K, std::void_t<decltype(std::declval<A &>()[std::declval<K>()])>>
        : std::true_type {};

}

/* EnumArray **************************************************************************************/

TEST(EnumArray, SubscriptsByKey) {
    EnumArray<fruit_t, int, FRUIT_SLOTS> counts{};
    counts[fruit_t::PEAR] = 7;
    EXPECT_EQ(7, counts[fruit_t::PEAR]);
    EXPECT_EQ(0, counts[fruit_t::APPLE]);
    EXPECT_EQ(FRUIT_SLOTS, counts.size());
}

TEST(EnumArray, AggregateInitializesInDeclarationOrder) {
    EnumArray<fruit_t, std::string, FRUIT_SLOTS> names = {{"apple", "pear", "plum"}};
    EXPECT_EQ("apple", names[fruit_t::APPLE]);
    EXPECT_EQ("pear", names[fruit_t::PEAR]);
    EXPECT_EQ("plum", names[fruit_t::PLUM]);
}

TEST(EnumArray, ConstSubscriptReads) {
    const EnumArray<fruit_t, int, FRUIT_SLOTS> counts = {{1, 2, 3}};
    EXPECT_EQ(2, counts[fruit_t::PEAR]);
}

TEST(EnumArray, IteratesValuesInKeyOrder) {
    EnumArray<fruit_t, int, FRUIT_SLOTS> counts = {{10, 20, 30}};
    std::vector<int> seen;
    for (int v : counts) seen.push_back(v);
    EXPECT_EQ((std::vector<int>{10, 20, 30}), seen);
}

// the whole point of the type: a subscript from another family must not compile
TEST(EnumArray, RejectsForeignKeyAtCompileTime) {
    using arr = EnumArray<fruit_t, int, FRUIT_SLOTS>;
    EXPECT_TRUE((has_subscript<arr, fruit_t>::value));
    EXPECT_FALSE((has_subscript<arr, veg_t>::value));
    EXPECT_FALSE((has_subscript<arr, int>::value));
}

/* EnumRange **************************************************************************************/

TEST(EnumRange, VisitsEveryEnumeratorInOrder) {
    std::vector<int> seen;
    for (fruit_t f : EnumRange<fruit_t, FRUIT_SLOTS>{}) seen.push_back(static_cast<int>(f));
    EXPECT_EQ((std::vector<int>{0, 1, 2}), seen);
}

TEST(EnumRange, PairsWithEnumArraySubscript) {
    EnumArray<fruit_t, std::string, FRUIT_SLOTS> names = {{"apple", "pear", "plum"}};
    std::string joined;
    for (fruit_t f : EnumRange<fruit_t, FRUIT_SLOTS>{}) joined += names[f] + ",";
    EXPECT_EQ("apple,pear,plum,", joined);
}
