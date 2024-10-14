#include "gtest/gtest.h"
#include "../src/dist.h"
#include "../src/globals.h"

namespace {

TEST(CigarScoreTest, BasicCigarScoreCalc) {

    std::vector<int> cigar = {PTR_MAT, PTR_INS, PTR_INS, PTR_MAT};
    EXPECT_EQ(7, calc_cig_swg_score(cigar, 3 /* sub */, 5 /* open */, 2 /* extend */));
}

} // namespace
