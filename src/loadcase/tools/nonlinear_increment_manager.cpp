/**
 * @file nonlinear_increment_manager.cpp
 * @brief Implements common scheduling and retry decisions for nonlinear steps.
 *
 * Both load control and arc-length control require bounded increments,
 * iteration-dependent adaptation, rollback-aware cutbacks and reliable
 * termination diagnostics. This module implements those shared decisions
 * without interpreting a load factor or evaluating nonlinear equilibrium.
 *
 * Controllers retain ownership of trial-state transactions, Newton solves,
 * endpoint handling and path-specific predictor or corrector calculations.
 *
 * @see NonlinearIncrementManager
 * @see LoadControl
 * @see ArcLengthControl
 */

#include "nonlinear_increment_manager.h"

#include "../../core/logging.h"

#include <algorithm>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * Validates and initializes the increment scheduler for a new nonlinear step.
 *
 * The configured increment bounds apply to the scalar step measure chosen by
 * the owning controller. Resetting restores the initial increment, clears
 * accepted-increment and consecutive-retry counts, and forgets prior errors.
 * It does not change the configuration values themselves.
 */
void NonlinearIncrementManager::reset() {
    // Validate settings shared by load control and arc-length control
    logging::error(maximum_increments > 0,
        "NonlinearIncrementManager requires maximum_increments > 0");
    logging::error(maximum_cutbacks > 0,
        "NonlinearIncrementManager requires maximum_cutbacks > 0");
    logging::error(initial_increment > Precision(0),
        "NonlinearIncrementManager requires initial_increment > 0");
    logging::error(minimum_increment > Precision(0),
        "NonlinearIncrementManager requires minimum_increment > 0");
    logging::error(maximum_increment >= minimum_increment,
        "NonlinearIncrementManager requires maximum_increment >= minimum_increment");
    logging::error(initial_increment >= minimum_increment &&
                   initial_increment <= maximum_increment,
        "NonlinearIncrementManager requires initial_increment within increment limits");
    logging::error(growth_factor > Precision(0),
        "NonlinearIncrementManager requires growth_factor > 0");
    logging::error(cutback_factor > Precision(0) && cutback_factor < Precision(1),
        "NonlinearIncrementManager requires cutback_factor between 0 and 1");
    logging::error(severe_cutback_factor > Precision(0) &&
                   severe_cutback_factor <= cutback_factor,
        "NonlinearIncrementManager requires severe_cutback_factor between 0 and cutback_factor");
    logging::error(fast_iterations > 0,
        "NonlinearIncrementManager requires fast_iterations > 0");
    logging::error(slow_iterations >= fast_iterations,
        "NonlinearIncrementManager requires slow_iterations >= fast_iterations");

    // Initialize the independent scheduling state
    increment_           = initial_increment;
    accepted_increments_ = 0;
    cutback_count_       = 0;
    failure_reason_      = "NONE";
}

Precision NonlinearIncrementManager::increment() const {
    return increment_;
}

Index NonlinearIncrementManager::accepted_increments() const {
    return accepted_increments_;
}

Index NonlinearIncrementManager::cutback_count() const {
    return cutback_count_;
}

bool NonlinearIncrementManager::can_attempt() const {
    return accepted_increments_ < maximum_increments;
}

const char* NonlinearIncrementManager::failure_reason() const {
    return failure_reason_.c_str();
}

/**
 * Replaces the planned increment without committing or rejecting a trial.
 *
 * The arc-length predictor may need to shorten its radius when approaching a
 * prescribed final load factor. This geometric correction does not itself
 * represent Newton convergence or failure and therefore does not affect counts.
 *
 * @param value New path increment size chosen by the owning controller.
 */
void NonlinearIncrementManager::set_increment(Precision value) {
    increment_ = value;
}

/**
 * Records one converged equilibrium increment and adapts the next step size.
 *
 * Only accepted increments consume maximum_increments. A successful solve
 * clears consecutive cutbacks, then scales the next step according to the
 * observed Newton iteration count. Fixed stepping restores the configured
 * initial increment, including after a shortened final-path proposal.
 *
 * @param newton_iterations Number of iterations in the converged Newton solve.
 */
void NonlinearIncrementManager::accept(Index newton_iterations) {
    // The owning controller has already committed its nonlinear trial state
    ++accepted_increments_;
    cutback_count_ = 0;

    // Fixed increments do not carry shortened endpoint attempts forward
    if (!adaptive) {
        increment_ = initial_increment;
        return;
    }

    // Enlarge quickly converged steps and reduce expensive converged steps
    if (newton_iterations <= fast_iterations) {
        increment_ *= growth_factor;
    } else if (newton_iterations >= slow_iterations) {
        increment_ *= cutback_factor;
    }

    // The new nominal increment must satisfy the configured size limits
    increment_ = std::clamp(increment_, minimum_increment, maximum_increment);
}

/**
 * Reduces the size of a failed nonlinear increment before another trial.
 *
 * Failed attempts are not counted as accepted increments. The reduction is
 * based on the increment actually attempted, not a difference of two rounded
 * load factors. Load control can use a stronger cutback for divergence, line
 * search failures or analysis exceptions; arc-length control uses the normal
 * factor. The owning controller must restore its solution and material/contact
 * trial state before invoking this method.
 *
 * @param attempted_increment Nominal or endpoint-adjusted attempted step.
 * @param severe Selects the stronger reduction factor when true.
 * @param detail Optional error description from the rejected solve.
 * @return True when another trial is allowed; false when stepping must stop.
 */
bool NonlinearIncrementManager::reject(Precision attempted_increment,
                                       bool      severe,
                                       const std::string& detail) {
    // Direct/fixed nonlinear control must fail instead of attempting a cutback
    if (!adaptive) {
        failure_reason_ = detail.empty()
            ? "FIXED_INCREMENT_FAILED"
            : "FIXED_INCREMENT_FAILED: " + detail;
        return false;
    }

    const Precision factor = severe ? severe_cutback_factor : cutback_factor;
    const bool can_retry = retry_with(attempted_increment * factor);

    // Retain the underlying exception text when a limit stopped the retry
    if (!can_retry && !detail.empty()) {
        failure_reason_ += ": " + detail;
    }

    return can_retry;
}

/**
 * Schedules a replacement increment following a rejected trial.
 *
 * This method also supports arc-length final-step adjustment, which can
 * request a retry despite an otherwise fixed increment policy. Every retry
 * consumes one cutback allowance and must respect the lower increment limit.
 *
 * @param value New increment size requested by the owning controller.
 * @param maximum_failure_reason Failure classification when retries are spent.
 * @return True if the requested retry is permitted.
 */
bool NonlinearIncrementManager::retry_with(Precision value, const char* maximum_failure_reason) {
    // Count controller-requested retries and install the proposed step
    increment_ = value;
    ++cutback_count_;

    // Do not silently continue below the configured minimum or retry budget
    if (increment_ < minimum_increment) {
        failure_reason_ = "MINIMUM_INCREMENT";
        return false;
    }

    if (cutback_count_ > maximum_cutbacks) {
        failure_reason_ = maximum_failure_reason;
        return false;
    }

    return true;
}

void NonlinearIncrementManager::fail(const std::string& reason) {
    failure_reason_ = reason;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
