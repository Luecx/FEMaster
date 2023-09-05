#pragma once

#include "types.h"
#include "logging.h"

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
    static auto measure(Func&& func, const std::string& description) {
        Timer timer;
        logging::info(true, "");
        logging::info(true, "Begin of ", std::setw(75), std::left, description);
        timer.start();
        logging::up();
        auto result = func();
        logging::down();
        timer.stop();
        logging::info(true, "Finished ",
                 std::setw(75), std::left, description, "[",
                 std::setw(6) , std::right, std::to_string(timer.elapsed()), " ms]");
        return result;
    }

    private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
};