/**
 * @file timer.cpp
 * @brief Implements the timer helper used for performance measurements.
 *
 * @see src/core/timer.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "timer.h"

namespace fem {

void Timer::start() {
    start_time = Clock::now();
}

void Timer::stop() {
    end_time = Clock::now();
}

Time Timer::elapsed() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
}
} // namespace fem
