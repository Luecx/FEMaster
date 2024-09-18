#pragma once

#include "logging.h"
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

    template<typename Func>
    static auto measure(Func&& func, const std::string& description, bool log_output = true) {
        Timer timer;

        // If logging is enabled, print the start message
        if (log_output) {
            logging::info(true, "");
            logging::info(true, "Begin of ", std::setw(75), std::left, description);
            logging::up();
        }

        // Start the timer
        timer.start();

        // Check if the function returns void
        if constexpr (std::is_void_v<decltype(func())>) {
            func();          // Execute the function if it returns void
            timer.stop();    // Stop the timer
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
            // Execute the function and return the result if not void
            auto result = func();
            timer.stop();    // Stop the timer

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

            return result;    // Return the result of the function
        }
    }

    template<typename Func>
    static double measure_time(Func&& func) {
        Timer timer;

        // Start the timer
        timer.start();

        // Function is assumed to return void, we measure time and return it
        func();    // Execute the function

        // Stop the timer and return the elapsed time
        timer.stop();
        return timer.elapsed();
    }

    private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
};
