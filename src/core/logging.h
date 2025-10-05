/**
 * @file logging.h
 * @brief Declares lightweight logging helpers used throughout the project.
 *
 * Provides indentation management and templated logging utilities implemented
 * in `logging.ipp`.
 *
 * @see src/core/logging.ipp
 * @see src/core/timer.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include <string>

namespace fem {
namespace logging {

extern int indentation_level; ///< Tracks the current indentation depth for log output.

/**
 * @brief Increases the indentation used for subsequent log messages.
 */
void up();

/**
 * @brief Decreases the indentation level if possible.
 */
void down();

/**
 * @brief Returns the indentation prefix string associated with the current level.
 */
std::string get_indentation();

/**
 * @brief Logs a warning when the supplied condition evaluates to false.
 */
template<typename... Args>
void warning(bool condition, Args... args);

/**
 * @brief Logs an informational message when the supplied condition is true.
 */
template<typename... Args>
void info(bool condition, Args... args);

/**
 * @brief Logs an error and aborts when the supplied condition evaluates to false.
 */
template<typename... Args>
void error(bool condition, Args... args);

} // namespace logging
} // namespace fem

#include "logging.ipp"
