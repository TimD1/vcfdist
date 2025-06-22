#include "gtest/gtest.h"
#include "../../../src/dist.h"
#include "../../../src/globals.h"

namespace {

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
