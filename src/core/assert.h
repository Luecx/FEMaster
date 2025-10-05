/**
 * @file assert.h
 * @brief Declares runtime assertion and check macros used across the codebase.
 *
 * The helpers provide lightweight diagnostics that raise exceptions when
 * conditions are violated. In release builds (`NDEBUG` defined) the assertions
 * remain active, whereas runtime checks always execute regardless of the build
 * mode.
 *
 * @see src/core/logging.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include <iostream>
#include <stdexcept>

#ifdef NDEBUG
/**
 * @brief Verifies a condition and throws if it evaluates to false (active in release builds).
 *
 * @param condition Boolean expression to evaluate.
 * @param message Diagnostic message emitted and used for the thrown exception.
 */
#define runtime_assert(condition, message)                                                                             \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cout << "runtime assertion triggered on: " << #condition << "\n";                                     \
            std::cout << "message: " << message << "\n";                                                               \
            std::cout << "line   : " << __LINE__ << "\n";                                                              \
            std::cout << "file   : " << __FILE__ << "\n";                                                              \
            throw std::runtime_error(message);                                                                         \
        }                                                                                                              \
    } while (false)
#else
#define runtime_assert(condition, message)
#endif

/**
 * @brief Always-on condition check that throws when the supplied expression fails.
 *
 * @param condition Boolean expression to evaluate.
 * @param message Diagnostic message emitted and used for the thrown exception.
 */
#define runtime_check(condition, message)                                                                              \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cout << "runtime check triggered on: " << #condition << "\n";                                         \
            std::cout << "message: " << message << "\n";                                                               \
            std::cout << "line   : " << __LINE__ << "\n";                                                              \
            std::cout << "file   : " << __FILE__ << "\n";                                                              \
            throw std::runtime_error(message);                                                                         \
        }                                                                                                              \
    } while (false)
