/******************************************************************************
* @file timer.h
* @brief Provides functionality for timing code execution in the FEMaster solver.
*
* This class allows for measuring the time taken to execute functions, either by
* logging the results or returning the elapsed time in milliseconds. It supports
* both void and non-void function calls and logs the execution time if desired.
*
* The `Timer` class can be used to measure time by calling `start()` and `stop()`,
* or by using the `measure` function to time a callable and optionally log the result.
*
* @note The class uses high-resolution clocks from the `<chrono>` library to ensure
* precise time measurements.
*
* @author Finn Eggers
* @date 24.10.2024
******************************************************************************/

#pragma once  // Ensures this file is only included once during compilation

#include "logging.h"  // For logging timing information
#include "types.h"    // For Time typedef

#include <chrono>     // For high-resolution timing

namespace fem {

class Timer {
   public:
   /******************************************************************************
    * @brief Starts the timer by capturing the current time.
    *
    * This function stores the current time using a high-resolution clock, marking
    * the beginning of a timed section.
    ******************************************************************************/
   void start();

   /******************************************************************************
    * @brief Stops the timer by capturing the current time.
    *
    * This function stores the current time, marking the end of the timed section.
    ******************************************************************************/
   void stop();

   /******************************************************************************
    * @brief Returns the elapsed time between the start and stop in milliseconds.
    *
    * This function calculates the time elapsed between the last call to `start()`
    * and `stop()`, returning it as a `Time` (defined as `uint64_t`).
    *
    * @return Time Elapsed time in milliseconds.
    ******************************************************************************/
   [[nodiscard]] Time elapsed() const;

   /******************************************************************************
    * @brief Measures the time taken to execute a function and optionally logs it.
    *
    * This template function times the execution of a provided function, logs the
    * description and time taken (in milliseconds) if logging is enabled, and
    * returns the result of the function.
    *
    * @tparam Func A callable type (e.g., lambda, function).
    * @param func The function to execute and time.
    * @param description A string describing the timed operation.
    * @param log_output A boolean flag indicating whether to log the output (default: true).
    * @return If `func` is non-void, returns the result of the function. If `func` is void, nothing is returned.
    ******************************************************************************/
   template<typename Func>
   static auto measure(Func&& func, const std::string& description, bool log_output = true) {
       Timer timer;

       // If logging is enabled, print the start message
       if (log_output) {
           logging::info(true, "");  // Log empty line for clarity
           logging::info(true, "Begin of ", std::setw(75), std::left, description);  // Log the start of the operation
           logging::up();  // Increase log indentation
       }

       // Start the timer
       timer.start();

       // Check if the function returns void
       if constexpr (std::is_void_v<decltype(func())>) {
           func();          // Execute the function if it returns void
           timer.stop();    // Stop the timer

           if (log_output) {
               logging::down();  // Decrease log indentation
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
                             " ms]");  // Log the elapsed time
           }
       } else {
           // Execute the function and return the result if it does not return void
           auto result = func();
           timer.stop();    // Stop the timer

           if (log_output) {
               logging::down();  // Decrease log indentation
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
                             " ms]");  // Log the elapsed time
           }

           return result;  // Return the result of the function
       }
   }

   /******************************************************************************
    * @brief Measures the time taken to execute a function and returns it.
    *
    * This function times the execution of a provided function (assumed to return void)
    * and returns the elapsed time in milliseconds.
    *
    * @tparam Func A callable type (e.g., lambda, function).
    * @param func The function to execute and time.
    * @return double The time taken by the function in milliseconds.
    ******************************************************************************/
   template<typename Func>
   static double measure_time(Func&& func) {
       Timer timer;

       // Start the timer
       timer.start();

       // Execute the function (assumed to return void)
       func();

       // Stop the timer and return the elapsed time
       timer.stop();
       return timer.elapsed();
   }

   private:
   std::chrono::high_resolution_clock::time_point start_time;  ///< The start time of the timer
   std::chrono::high_resolution_clock::time_point end_time;    ///< The end time of the timer
};

}  // namespace fem