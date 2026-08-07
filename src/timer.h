/**
 * @file timer.h
 * @brief Wall-clock timer for measuring elapsed time of named pipeline stages.
 */
#ifndef _TIMER_H_
#define _TIMER_H_

#include <string>
#include <chrono>

#include "defs.h"

/**
 * @class timer
 * @brief Measures and accumulates elapsed wall-clock time for a named pipeline stage.
 */
class timer {
public:
    /** @brief Constructs timer with default name "default". */
    timer() {};

    /** @brief Constructs timer with specified name. */
    timer(const std::string & timer_name) : name(timer_name) {};

    /** @brief Records current system time and marks timer as running. */
    void start();

    /** @brief Returns the timer's name string. */
    std::string get_name();

    /** @brief Stops timer and accumulates elapsed time since last start(). */
    void stop();

    /** @brief Returns total accumulated elapsed time in seconds. */
    double total();

    /** @brief Prints timer index, name, and total elapsed time to console. */
    void print(int i);

private:
    std::string name = "default"; ///< Timer label string
    std::chrono::time_point<std::chrono::system_clock> start_time; ///< System clock time point recorded at last start()
    double total_time = 0;        ///< Accumulated elapsed time in nanoseconds
    bool running = false;         ///< True if timer is currently running
};

/** @brief Writes all pipeline stage timer names and elapsed times to TSV file. */
void write_runtime();

#endif
