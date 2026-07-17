/**
 * @file load_control.cpp
 * @brief Implements a generic nonlinear load-control solver.
 */

#include "load_control.h"

#include "../../core/logging.h"

#include <algorithm>

namespace fem {
namespace loadcase {
namespace tools {

bool LoadControl::solve(
    DynamicVector&           q,
    Precision&               lambda,
    const Evaluate&          evaluate,
    const LinearSolve&       linear_solve,
    const ResidualNorm&      residual_norm,
    const CorrectionNorm&    correction_norm,
    const IterationCallback& on_iteration,
    const IncrementCallback& on_increment
) {
    logging::error(maximum_increments > 0,
        "LoadControl requires maximum_increments > 0");
    logging::error(maximum_iterations > 0,
        "LoadControl requires maximum_iterations > 0");
    logging::error(tolerance > Precision(0),
        "LoadControl requires tolerance > 0");
    logging::error(initial_increment > Precision(0),
        "LoadControl requires initial_increment > 0");
    logging::error(minimum_increment > Precision(0),
        "LoadControl requires minimum_increment > 0");
    logging::error(maximum_increment >= minimum_increment,
        "LoadControl requires maximum_increment >= minimum_increment");
    logging::error(initial_increment >= minimum_increment &&
                   initial_increment <= maximum_increment,
        "LoadControl requires initial_increment between minimum_increment and maximum_increment");
    logging::error(growth_factor > Precision(0),
        "LoadControl requires growth_factor > 0");
    logging::error(cutback_factor > Precision(0) && cutback_factor < Precision(1),
        "LoadControl requires cutback_factor between 0 and 1");
    logging::error(fast_iterations > 0,
        "LoadControl requires fast_iterations > 0");
    logging::error(slow_iterations >= fast_iterations,
        "LoadControl requires slow_iterations >= fast_iterations");
    logging::error(maximum_cutbacks > 0,
        "LoadControl requires maximum_cutbacks > 0");

    logging::error(static_cast<bool>(evaluate),
        "LoadControl requires an evaluate callback");
    logging::error(static_cast<bool>(linear_solve),
        "LoadControl requires a linear solve callback");
    logging::error(static_cast<bool>(residual_norm),
        "LoadControl requires a residual norm callback");
    logging::error(static_cast<bool>(correction_norm),
        "LoadControl requires a correction norm callback");

    reset_state_();

    Index cutback_count = 0;

    while (lambda < Precision(1) - tolerance &&
           accepted_increments_ < maximum_increments) {
        increment_ = std::min(increment_, Precision(1) - lambda);

        if (increment_ < minimum_increment) {
            failure_reason_ = "MINIMUM_INCREMENT";
            return false;
        }

        const DynamicVector q_accepted      = q;
        const Precision     lambda_accepted = lambda;
        const Precision     target_lambda   = lambda + increment_;

        configure_newton_();

        const bool converged = newton_.solve(
            q,
            [&](const DynamicVector& current_q,
                DynamicVector&       residual,
                SparseMatrix&        tangent) {
                evaluate(current_q, target_lambda, residual, tangent);
            },
            linear_solve,
            [&](const DynamicVector& residual) {
                return residual_norm(residual, target_lambda);
            },
            correction_norm,
            [&](Index     iteration,
                Precision current_residual_norm,
                Precision current_correction_norm,
                Precision convergence_order,
                Time      assembly_ms,
                Time      solve_ms,
                bool      iteration_converged) {
                if (on_iteration) {
                    on_iteration(
                        accepted_increments_ + 1,
                        iteration,
                        target_lambda,
                        current_residual_norm,
                        current_correction_norm,
                        convergence_order,
                        assembly_ms,
                        solve_ms,
                        iteration_converged
                    );
                }
            }
        );

        if (!converged) {
            q      = q_accepted;
            lambda = lambda_accepted;

            if (!adaptive) {
                failure_reason_ = "FIXED_INCREMENT_FAILED";
                return false;
            }

            increment_ *= cutback_factor;
            ++cutback_count;

            logging::info(true,
                "Increment rejected at lambda = ",
                target_lambda,
                "; reducing increment to ",
                increment_
            );

            if (increment_ < minimum_increment) {
                failure_reason_ = "MINIMUM_INCREMENT";
                return false;
            }

            if (cutback_count > maximum_cutbacks) {
                failure_reason_ = "MAXIMUM_CUTBACKS";
                return false;
            }

            continue;
        }

        lambda = target_lambda;

        ++accepted_increments_;
        cutback_count = 0;

        if (on_increment) {
            on_increment(accepted_increments_, q, lambda);
        }

        adapt_increment_();

        logging::info(true,
            "Accepted increment "   , accepted_increments_,
            ": lambda = "           , lambda,
            ", Newton iterations = ", newton_.iterations(),
            ", next increment = "  , increment_
        );
    }

    if (lambda < Precision(1) - tolerance) {
        failure_reason_ = "MAXIMUM_INCREMENTS";
        return false;
    }

    return true;
}

Index LoadControl::accepted_increments() const {
    return accepted_increments_;
}

Precision LoadControl::increment() const {
    return increment_;
}

const char* LoadControl::failure_reason() const {
    return failure_reason_;
}

void LoadControl::reset_state_() {
    accepted_increments_ = 0;
    increment_           = initial_increment;
    failure_reason_      = "NONE";
}

void LoadControl::configure_newton_() {
    newton_.maximum_iterations       = maximum_iterations;
    newton_.residual_tolerance       = tolerance;
    newton_.stagnation_tolerance     = Precision(1e-3) * tolerance;
    newton_.check_finite             = true;
    newton_.early_failure_detection  = true;
}

void LoadControl::adapt_increment_() {
    if (!adaptive) {
        increment_ = initial_increment;
        return;
    }

    if (newton_.iterations() <= fast_iterations) {
        increment_ *= growth_factor;
    } else if (newton_.iterations() >= slow_iterations) {
        increment_ *= cutback_factor;
    }

    increment_ = std::clamp(
        increment_,
        minimum_increment,
        maximum_increment
    );
}

} // namespace tools
} // namespace loadcase
} // namespace fem
