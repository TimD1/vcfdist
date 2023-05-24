#ifndef _TIMER_H_
#define _TIMER_H_

#include <string>
#include <chrono>

#include "globals.h"
#include "defs.h"

class timer {
public:
    timer() {};
    timer(std::string timer_name) : name(timer_name) {};

    void start();
    void stop();
    double total();
    void print(int depth);

private:
    std::string name = "default";
    std::chrono::time_point<std::chrono::system_clock> start_time;
    double total_time = 0;
    bool running = false;
};

#endif
