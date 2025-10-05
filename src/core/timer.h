/******************************************************************************
 * @file timer.h
 * @brief Declares utilities for timing code execution.
 *
 * Provides a small `Timer` helper with RAII-style measurement helpers that
 * integrate with the logging subsystem.
 *
 * @see src/core/timer.cpp
 * @see src/core/logging.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "logging.h"
#include "types_eig.h"

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>

namespace fem {

/******************************************************************************
 * @class Timer
 * @brief Measures elapsed wall-clock time.
 ******************************************************************************/
class Timer {
public:
    /// Captures the current time and starts the timer.
    void start();

    /// Captures the current time and stops the timer.
    void stop();

    /// Returns the elapsed time in milliseconds between `start` and `stop`.
    [[nodiscard]] Time elapsed() const;

    /******************************************************************************
     * @brief Measures a callable and optionally logs the duration.
     *
     * @tparam Func Callable type.
     * @param func Callable to invoke.
     * @param description Textual description used for logging.
     * @param log_output Whether to emit timing information via logging.
     * @return The callable's return value if non-void.
     ******************************************************************************/
    template<typename Func>
    static auto measure(Func&& func, const std::string& description, bool log_output = true);

    /******************************************************************************
     * @brief Measures a callable and returns only the elapsed time.
     *
     * @tparam Func Callable type.
     * @param func Callable to invoke.
     * @return Elapsed time in milliseconds.
     ******************************************************************************/
    template<typename Func>
    static double measure_time(Func&& func);

private:
    std::chrono::high_resolution_clock::time_point start_time{};
    std::chrono::high_resolution_clock::time_point end_time{};
};

template<typename Func>
auto Timer::measure(Func&& func, const std::string& description, bool log_output) {
    Timer timer;
    if (log_output) {
        logging::info(true, "");
        logging::info(true, "Begin of ", std::setw(75), std::left, description);
        logging::up();
    }

    timer.start();

    if constexpr (std::is_void_v<std::invoke_result_t<Func&>>) {
        std::forward<Func>(func)();
        timer.stop();
        if (log_output) {
            logging::down();
            logging::info(true,
                          "Finished ",
                          std::setw(75),
                          std::left,
                          description,
                          "[",
                          std::setw(6),
                          std::right,
                          std::setprecision(4),
                          timer.elapsed(),
                          " ms]");
        }
    } else {
        auto result = std::forward<Func>(func)();
        timer.stop();
        if (log_output) {
            logging::down();
            logging::info(true,
                          "Finished ",
                          std::setw(75),
                          std::left,
                          description,
                          "[",
                          std::setw(6),
                          std::right,
                          std::setprecision(4),
                          timer.elapsed(),
                          " ms]");
        }
        return result;
    }
}

template<typename Func>
double Timer::measure_time(Func&& func) {
    Timer timer;
    timer.start();
    func();
    timer.stop();
    return timer.elapsed();
}

} // namespace fem
