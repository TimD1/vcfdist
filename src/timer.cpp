#include "timer.h"
#include "defs.h"

void timer::start() {
    if (running) {
        ERROR("Cannot start an already-running timer.");
    }
    start_time = std::chrono::system_clock::now();
    running = true;
}

void timer::stop() {
    if (!running) {
        ERROR("Cannot stop an already-stopped timer.");
    }
    auto stop_time = std::chrono::system_clock::now();
    total_time += (std::chrono::duration_cast<std::chrono::milliseconds>(
            stop_time-start_time).count()) / 1000.0;
    running = false;
}

double timer::total() { 
    if (running) {
        ERROR("Must stop timer before calculating total time.");
    }
    return total_time; 
}

void timer::print() {
    INFO("%15s timer: %8.3fs", name.data(), total());
}
