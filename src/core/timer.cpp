#include "timer.h"

void Timer::start() {
    start_time = std::chrono::high_resolution_clock::now();
}

void Timer::stop() {
    end_time = std::chrono::high_resolution_clock::now();
}

uint64_t Timer::elapsed() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
}