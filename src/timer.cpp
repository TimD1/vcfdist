#include "timer.h"
#include "defs.h"

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
    INFO("  [%d] %-16s: %8.3fs", i, name.data(), total());
}
