/******************************************************************************
* @file timer.cpp
* @brief Implements the Timer class functions for measuring code execution time.
*
* This source file contains the implementation of the Timer class functions,
* which start and stop a high-resolution clock and calculate the elapsed time
* in milliseconds between the start and stop times.
*
* @note The Timer class uses the high-resolution clock from the `<chrono>` library
* for accurate time measurements.
*
* @see timer.h
* @author Finn Eggers
* @date 24.10.2024
******************************************************************************/

#include "timer.h"

/******************************************************************************
* @brief Starts the timer by capturing the current high-resolution clock time.
*
* This function stores the current time from the high-resolution clock in the
* `start_time` member variable, marking the beginning of a timed section.
******************************************************************************/
void fem::Timer::start() {
   start_time = std::chrono::high_resolution_clock::now();  // Capture the current time
}

/******************************************************************************
* @brief Stops the timer by capturing the current high-resolution clock time.
*
* This function stores the current time from the high-resolution clock in the
* `end_time` member variable, marking the end of a timed section.
******************************************************************************/
void fem::Timer::stop() {
   end_time = std::chrono::high_resolution_clock::now();  // Capture the current time
}

/******************************************************************************
* @brief Returns the elapsed time in milliseconds between the start and stop times.
*
* This function calculates the duration between the `start_time` and `end_time`,
* and returns the result in milliseconds. If the timer was never stopped after
* starting, the elapsed time may not be accurate.
*
* @return uint64_t The elapsed time in milliseconds.
******************************************************************************/
uint64_t fem::Timer::elapsed() const {
   // Calculate the difference between start and end times, and convert to milliseconds
   return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
}
