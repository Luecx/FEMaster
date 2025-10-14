/**
 * @file logging.h
 * @brief Declares lightweight logging helpers used throughout the project.
 *
 * Provides indentation management, maximum console width, and templated logging
 * utilities implemented in `logging.ipp`.
 *
 * Wrapping rules:
 *  - Hard newlines '\n' are respected.
 *  - Lines are wrapped to the configured console width.
 *  - On wrapped lines, the same prefix ([INFO]/[WARNING]/[ERROR] + indentation)
 *    is repeated.
 *  - Leading spaces/tabs of the original text segment (after the prefix) are
 *    detected and re-applied to all wrapped continuation lines to keep the
 *    visual left alignment of the text block.
 *
 * @see src/core/logging.ipp
 * @see src/core/timer.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include <cstddef>
#include <string>

namespace fem {
namespace logging {

extern int indentation_level; ///< Tracks the current indentation depth for log output.

/**
 * @brief Set the assumed console width (in characters) for wrapping.
 *        Default is 100. Must be >= 20 to avoid pathological behavior.
 */
void set_console_width(std::size_t cols);

/**
 * @brief Get the currently configured console width.
 */
std::size_t get_console_width();

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
