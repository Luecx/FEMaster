/**
 * @file run_quiet.h
 * @brief Provides callback execution with temporarily disabled log output.
 *
 * This core utility suppresses informational and warning output while a caller
 * executes a callback. Runtime error checks remain active in the logging
 * subsystem. Output that was enabled on entry is re-enabled on both normal
 * return and exception propagation; the caller owns the callback's work.
 *
 * @see fem::logging
 * @author Finn Eggers
 * @date 10.11.2026
 */

#pragma once

#include "logging.h"

#include <type_traits>
#include <utility>

namespace fem {

/**
 * Executes a callback without informational or warning logging.
 *
 * Save the logging state before invoking the perfectly forwarded callable.
 * Re-enable output afterward if it was enabled on entry, including when the
 * callable throws. Void callbacks return no value; other results are returned
 * by value, matching the existing solver wrapper's behavior.
 *
 * The callback should leave logging disabled while it runs. Logging state is
 * global, so suppression is shared with other work using the same logger.
 *
 * @param function Callback executed once while output is disabled.
 * @return Callback result by value, or void for a void callback.
 */
template<class Function>
decltype(auto) run_quiet(Function&& function) {
    // Suppress output without disabling runtime error checks.
    const bool logging_was_enabled = logging::is_enabled();
    logging::disable();

    // Restore previously enabled output on every return or exception path.
    try {
        if constexpr (std::is_void_v<std::invoke_result_t<Function>>) {
            std::forward<Function>(function)();
            if (logging_was_enabled) logging::enable();
            return;
        } else {
            auto result = std::forward<Function>(function)();
            if (logging_was_enabled) logging::enable();
            return result;
        }
    } catch (...) {
        if (logging_was_enabled) logging::enable();
        throw;
    }
}

} // namespace fem
