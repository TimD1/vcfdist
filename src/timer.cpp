#include "timer.h"
#include "defs.h"
#include "globals.h"

void timer::start() {
    if (running) {
        ERROR("Cannot start an already-running timer (%s).", name.data());
    }
    start_time = std::chrono::system_clock::now();
    running = true;
}

void timer::stop() {
    if (!running) {
        ERROR("Cannot stop an already-stopped timer (%s).", name.data());
    }
    auto stop_time = std::chrono::system_clock::now();
    total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(
            stop_time-start_time).count();
    running = false;
}

double timer::total() { 
    if (running) {
        ERROR("Must stop timer before calculating total time (%s).", name.data());
    }
    return total_time / 1000000000.0; 
}

void timer::print(int i) {
    INFO("  [%d] %-37s: %8.3fs", i, name.data(), total());
}

std::string timer::get_name() {
    return name;
}

void write_runtime() {
    std::string runtimes_fn = g.out_prefix + "runtime.tsv";
    if (g.verbosity >= 1) INFO("  Writing stage runtimes '%s'", runtimes_fn.data());
    FILE* out_runtimes = fopen(runtimes_fn.data(), "w");
    for (int i = 0; i <= TIME_TOTAL; i++) {
        fprintf(out_runtimes, "%s\t%lf\n", g.timers[i].get_name().data(), g.timers[i].total());
    }
    fclose(out_runtimes);
}
