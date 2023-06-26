#pragma once

#include "types.h"

#include <chrono>

class Timer {
    public:
    // starts the timer
    void start();

    // stops the timer
    void stop();

    // returns the elapsed time in milliseconds
    [[nodiscard]] Time elapsed() const;

    private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
};