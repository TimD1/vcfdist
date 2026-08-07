/**
 * @file timer.cpp
 * @brief Wall-clock timer implementation for named pipeline stages.
 */
#include "timer.h"
#include "defs.h"
#include "globals.h"

/**
 * @brief Records current system time and marks timer as running.
 * @throws Error if timer is already running.
 */
void timer::start() {
    if (running) {
        ERROR("Cannot start an already-running timer (%s).", name.data());
    }
    start_time = std::chrono::system_clock::now();
    running = true;
}

/**
 * @brief Stops timer and accumulates elapsed time since last start().
 * @throws Error if timer is not running.
 */
void timer::stop() {
    if (!running) {
        ERROR("Cannot stop an already-stopped timer (%s).", name.data());
    }
    auto stop_time = std::chrono::system_clock::now();
    total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(
            stop_time-start_time).count();
    running = false;
}

/**
 * @brief Returns total accumulated elapsed time in seconds.
 * @return Total elapsed seconds across all start/stop intervals
 * @throws Error if timer is still running.
 */
double timer::total() {
    if (running) {
        ERROR("Must stop timer before calculating total time (%s).", name.data());
    }
    return total_time / 1000000000.0;
}

/**
 * @brief Prints timer index, name, and total elapsed time to console.
 * @param[in] i Timer index used for labeling output
 */
void timer::print(int i) {
    INFO("  [%d] %-17s: %8.3fs", i, name.data(), total());
}

/**
 * @brief Returns the timer's name string.
 * @return Timer name provided at construction
 */
std::string timer::get_name() {
    return name;
}

/**
 * @brief Writes all pipeline stage timer names and elapsed times to TSV file.
 * @throws ERROR if the output runtime TSV file cannot be opened for writing
 */
void write_runtime() {
    std::string runtimes_fn = g.out_prefix + "runtime.tsv";
    if (g.verbosity >= 1) INFO("  Writing stage runtimes '%s'", runtimes_fn.data());
    FILE* out_runtimes = fopen(runtimes_fn.data(), "w");
    if (out_runtimes == NULL) {
        ERROR("Failed to open runtime TSV file '%s'", runtimes_fn.data());
    }
    for (timer_t t : EnumRange<timer_t, TIMER_SLOTS>{}) {
        fprintf(out_runtimes, "%s\t%lf\n", g.stage(t).get_name().data(), g.stage(t).total());
    }
    fclose(out_runtimes);
}
