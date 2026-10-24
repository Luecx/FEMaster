/**
 * @file nonlinear_increment_manager.h
 * @brief Declares shared increment scheduling for nonlinear solution procedures.
 *
 * NonlinearIncrementManager owns the increment size, adaptive growth and
 * cutback decisions, accepted-increment counts and termination diagnostics
 * shared by load control and arc-length control.
 *
 * The increment is a method-dependent scalar: it denotes a load-factor
 * increment in LoadControl and a normalized path-radius increment in
 * ArcLengthControl. The manager does not evaluate equilibrium equations,
 * choose path directions, calculate load factors, run Newton iterations or
 * own material and contact trial states. Those responsibilities remain with
 * the nonlinear controllers and NonlinearStatic.
 *
 * @see LoadControl
 * @see ArcLengthControl
 * @see NewtonSolver
 */

#pragma once

#include "../../core/types_num.h"

#include <string>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * @brief Manages attempted and accepted increments of a nonlinear analysis.
 *
 * The manager provides one common scheduling policy for load-controlled and
 * arc-length-controlled equilibrium paths. Its increment is a positive scalar
 * with prescribed initial, minimum and maximum values. Individual controllers
 * retain the interpretation of this scalar and decide whether a nonlinear
 * equilibrium point constitutes successful progress along their path.
 *
 * Only successfully converged increments increase accepted_increments_.
 * Failed attempts may reduce the next increment and increase cutback_count_.
 * The consecutive cutback count is reset following a successful increment,
 * while adaptive growth or reduction uses the accepted Newton iteration count.
 *
 * The manager deliberately does not store the finite-element solution or
 * commit nonlinear material/contact history. Controllers perform trial state
 * transactions and call accept, reject or retry_with as appropriate.
 */
class NonlinearIncrementManager {
public:
    // Limits on accepted increments and consecutive retries
    Index maximum_increments = 100;
    Index maximum_cutbacks   = 20;

    // Allowed initial and adaptive increment sizes. LoadControl interprets
    // these values as load-factor changes, ArcLengthControl as path radii.
    Precision initial_increment = Precision(0.1);
    Precision minimum_increment = Precision(1e-4);
    Precision maximum_increment = Precision(1);

    // Adaptive step control from the last accepted Newton iteration count
    Precision growth_factor         = Precision(1.5);
    Precision cutback_factor        = Precision(0.5);
    Precision severe_cutback_factor = Precision(0.25);

    Index fast_iterations = 6;
    Index slow_iterations = 10;

    bool adaptive = true;

private:
    // Current scheduling state, independent of the nonlinear solution vector
    Precision increment_ = Precision(0);

    Index accepted_increments_ = 0;
    Index cutback_count_       = 0;

    std::string failure_reason_ = "NONE";

public:
    // Construction
    NonlinearIncrementManager() = default;

    // Initialize and validate the increment policy for a fresh nonlinear step
    void reset();

    // Current increment and accumulated scheduling statistics
    Precision increment()           const;
    Index     accepted_increments() const;
    Index     cutback_count()       const;
    bool      can_attempt()         const;
    const char* failure_reason()    const;

    // Allow a controller to shorten a planned mathematical path increment.
    // In particular, arc-length control can trim its radius at the endpoint
    // without counting this change as a rejected nonlinear solve.
    void set_increment(Precision value);

    // Commit a converged increment and adapt the following increment from
    // the Newton iteration count. Only accepted increments consume the limit.
    void accept(Index newton_iterations);

    // Reject a failed nonlinear solve and prepare its next adaptive retry.
    // Severe failures may use a stronger cutback than ordinary Newton failure.
    bool reject(Precision attempted_increment, bool severe, const std::string& detail = {});

    // Reschedule a controller-specific retry such as arc-length endpoint
    // correction. Unlike reject, this may also be used with fixed increments.
    bool retry_with(Precision value, const char* maximum_failure_reason = "MAXIMUM_CUTBACKS");

    // Preserve a concrete termination reason for the owning controller
    void fail(const std::string& reason);
};

} // namespace tools
} // namespace loadcase
} // namespace fem
