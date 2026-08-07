/**
 * @file test_timer.cpp
 * @brief Unit tests for timer.cpp: named wall-clock timers and stage runtime output.
 */
#include <chrono>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/globals.h"
#include "../../../src/timer.h"
#include "test_helpers.h"

namespace {

/* timer ******************************************************************************************/

TEST(Timer, GetName) {
    GlobalsGuard guard;
    timer t("clustering");
    EXPECT_EQ("clustering", t.get_name());
}

TEST(Timer, DefaultName) {
    GlobalsGuard guard;
    timer t;
    EXPECT_EQ("default", t.get_name());
}

TEST(Timer, StartTwiceErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                timer t("clustering");
                t.start();
                t.start();
            }, testing::ExitedWithCode(1), "Cannot start an already-running timer");
}

TEST(Timer, StopNotRunningErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                timer t("clustering");
                t.stop();
            }, testing::ExitedWithCode(1), "Cannot stop an already-stopped timer");
}

TEST(Timer, TotalWhileRunningErrors) {
    GlobalsGuard guard;
    EXPECT_EXIT({
                timer t("clustering");
                t.start();
                t.total();
            }, testing::ExitedWithCode(1), "Must stop timer before calculating total time");
}

TEST(Timer, TotalAccumulates) {
    GlobalsGuard guard;
    timer t("clustering");

    // an unused timer has accumulated nothing
    EXPECT_EQ(0.0, t.total());

    t.start();
    t.stop();
    double one_interval = t.total();
    EXPECT_GE(one_interval, 0.0);

    // wall-clock durations are not reproducible, so only require monotonic accumulation
    t.start();
    t.stop();
    EXPECT_GE(t.total(), one_interval);
}

TEST(Timer, TotalSecondsScale) {
    GlobalsGuard guard;
    timer t("clustering");
    t.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    t.stop();

    // elapsed nanoseconds are divided by 1e9, so a 20ms interval is a small fraction of a second
    EXPECT_GE(t.total(), 0.005);
    EXPECT_LT(t.total(), 5.0);
}

TEST(Timer, Print) {
    GlobalsGuard guard;
    timer t("clustering");
    t.start();
    t.stop();

    testing::internal::CaptureStderr();
    t.print(idx(TIME_CLUSTER));
    std::string out = testing::internal::GetCapturedStderr();

    EXPECT_NE(std::string::npos, out.find("[" + std::to_string(idx(TIME_CLUSTER)) + "]"));
    EXPECT_NE(std::string::npos, out.find("clustering"));

    // the elapsed time is printed as "%8.3fs", so the line ends with a seconds suffix
    ASSERT_GE(out.size(), size_t(2));
    EXPECT_EQ('\n', out.back());
    EXPECT_EQ('s', out[out.size()-2]);
    EXPECT_NE(std::string::npos, out.find("."));
}

/* write_runtime **********************************************************************************/

TEST(WriteRuntime, OneRowPerStage) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";
    g.init_timers();

    write_runtime();

    std::ifstream in(g.out_prefix + "runtime.tsv");
    ASSERT_TRUE(in.is_open());
    std::vector<std::string> lines;
    std::string line;
    while (getline(in, line)) lines.push_back(line);

    ASSERT_EQ(idx(TIME_TOTAL)+1, lines.size());
    for (timer_t t : EnumRange<timer_t, TIMER_SLOTS>{}) {
        EXPECT_EQ(size_t(0), lines[idx(t)].rfind(timer_strs[t] + "\t", 0))
                << "row " << idx(t) << ": " << lines[idx(t)];
    }
}

TEST(WriteRuntime, UnwritableDirectoryErrors) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path("missing/");
    g.init_timers();

    EXPECT_EXIT(write_runtime(), testing::ExitedWithCode(1), "Failed to open runtime TSV file");
}

} // namespace
