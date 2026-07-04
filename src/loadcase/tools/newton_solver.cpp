/**
 * @file newton_solver.cpp
 * @brief Implements a generic Newton solver.
 */

#include "newton_solver.h"

#include "../../core/logging.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {
namespace loadcase {
namespace tools {

bool NewtonSolver::solve(
    DynamicVector&           x,
    const Evaluate&          evaluate,
    const LinearSolve&       linear_solve,
    const Norm&              residual_norm,
    const CorrectionNorm&    correction_norm,
    const IterationCallback& on_iteration
) {
    logging::error(maximum_iterations > 0,
        "NewtonSolver requires maximum_iterations > 0");
    logging::error(residual_tolerance > Precision(0),
        "NewtonSolver requires residual_tolerance > 0");
    logging::error(stagnation_tolerance >= Precision(0),
        "NewtonSolver requires stagnation_tolerance >= 0");
    logging::error(convergence_check_start > 0,
        "NewtonSolver requires convergence_check_start > 0");
    logging::error(divergence_factor > Precision(1),
        "NewtonSolver requires divergence_factor > 1");
    logging::error(minimum_residual_reduction >= Precision(0),
        "NewtonSolver requires minimum_residual_reduction >= 0");
    logging::error(stagnation_ratio >= Precision(0),
        "NewtonSolver requires stagnation_ratio >= 0");

    logging::error(static_cast<bool>(evaluate),
        "NewtonSolver requires an evaluate callback");
    logging::error(static_cast<bool>(linear_solve),
        "NewtonSolver requires a linear solve callback");
    logging::error(static_cast<bool>(residual_norm),
        "NewtonSolver requires a residual norm callback");
    logging::error(static_cast<bool>(correction_norm),
        "NewtonSolver requires a correction norm callback");

    reset_state_();

    DynamicVector residual;
    SparseMatrix  tangent;

    for (Index iter = 1; iter <= maximum_iterations; ++iter) {
        Timer assembly_timer;
        Timer solve_timer;

        assembly_timer.start();
        evaluate(x, residual, tangent);
        assembly_timer.stop();

        iterations_ = iter;

        update_residual_history_();

        last_residual_norm_ = residual_norm(residual);

        logging::error(
            !check_finite || std::isfinite(last_residual_norm_),
            "NewtonSolver: residual norm is NaN/Inf"
        );

        if (iter == 1) {
            initial_residual_norm_ = last_residual_norm_;
        }

        update_convergence_order_();
        update_failure_counters_();

        if (last_residual_norm_ <= residual_tolerance) {
            last_correction_norm_ = Precision(0);

            if (on_iteration) {
                on_iteration(
                    iter,
                    last_residual_norm_,
                    last_correction_norm_,
                    convergence_order_,
                    assembly_timer.elapsed(),
                    Time(0),
                    true
                );
            }

            return true;
        }

        if (should_stop_early_()) {
            if (on_iteration) {
                on_iteration(
                    iter,
                    last_residual_norm_,
                    last_correction_norm_,
                    convergence_order_,
                    assembly_timer.elapsed(),
                    Time(0),
                    false
                );
            }

            return false;
        }

        solve_timer.start();
        DynamicVector dx = linear_solve(tangent, residual);
        solve_timer.stop();

        logging::error(
            !check_finite || dx.allFinite(),
            "NewtonSolver: correction contains NaN/Inf entries"
        );

        last_correction_norm_ = correction_norm(x, dx);

        logging::error(
            !check_finite || std::isfinite(last_correction_norm_),
            "NewtonSolver: correction norm is NaN/Inf"
        );

        if (on_iteration) {
            on_iteration(
                iter,
                last_residual_norm_,
                last_correction_norm_,
                convergence_order_,
                assembly_timer.elapsed(),
                solve_timer.elapsed(),
                false
            );
        }

        x += dx;

        if (early_failure_detection &&
            iter >= convergence_check_start &&
            stagnation_tolerance > Precision(0) &&
            last_correction_norm_ <= stagnation_tolerance) {
            failed_by_stagnation_ = true;
            return false;
        }
    }

    failed_by_maximum_iterations_ = true;
    return false;
}

Index NewtonSolver::iterations() const {
    return iterations_;
}

Precision NewtonSolver::initial_residual_norm() const {
    return initial_residual_norm_;
}

Precision NewtonSolver::last_residual_norm() const {
    return last_residual_norm_;
}

Precision NewtonSolver::last_correction_norm() const {
    return last_correction_norm_;
}

Precision NewtonSolver::convergence_order() const {
    return convergence_order_;
}

bool NewtonSolver::failed_by_divergence() const {
    return failed_by_divergence_;
}

bool NewtonSolver::failed_by_residual_increase() const {
    return failed_by_residual_increase_;
}

bool NewtonSolver::failed_by_stagnation() const {
    return failed_by_stagnation_;
}

bool NewtonSolver::failed_by_poor_reduction() const {
    return failed_by_poor_reduction_;
}

bool NewtonSolver::failed_by_maximum_iterations() const {
    return failed_by_maximum_iterations_;
}

const char* NewtonSolver::failure_reason() const {
    if (failed_by_divergence_) {
        return "DIVERGENCE";
    }

    if (failed_by_residual_increase_) {
        return "RESIDUAL_INCREASE";
    }

    if (failed_by_stagnation_) {
        return "STAGNATION";
    }

    if (failed_by_poor_reduction_) {
        return "POOR_RESIDUAL_REDUCTION";
    }

    if (failed_by_maximum_iterations_) {
        return "MAXIMUM_ITERATIONS";
    }

    return "NONE";
}

void NewtonSolver::reset_state_() {
    iterations_ = 0;

    initial_residual_norm_           = Precision(0);
    last_residual_norm_              = Precision(0);
    last_correction_norm_            = Precision(0);
    previous_residual_norm_          = Precision(0);
    previous_previous_residual_norm_ = Precision(0);
    convergence_order_               = Precision(0);

    residual_increase_count_         = 0;
    stagnation_count_                = 0;

    failed_by_divergence_            = false;
    failed_by_residual_increase_     = false;
    failed_by_stagnation_            = false;
    failed_by_poor_reduction_        = false;
    failed_by_maximum_iterations_    = false;
}

void NewtonSolver::update_residual_history_() {
    previous_previous_residual_norm_ = previous_residual_norm_;
    previous_residual_norm_          = last_residual_norm_;
}

void NewtonSolver::update_convergence_order_() {
    convergence_order_ = Precision(0);

    if (iterations_ < 3) {
        return;
    }

    const Precision r0 = previous_previous_residual_norm_;
    const Precision r1 = previous_residual_norm_;
    const Precision r2 = last_residual_norm_;

    if (r0 <= Precision(0) || r1 <= Precision(0) || r2 <= Precision(0)) {
        return;
    }

    const Precision numerator   = std::log(r2 / r1);
    const Precision denominator = std::log(r1 / r0);

    if (std::abs(denominator) <= Precision(1e-14)) {
        return;
    }

    convergence_order_ = numerator / denominator;
}

void NewtonSolver::update_failure_counters_() {
    if (iterations_ <= 1) {
        return;
    }

    if (previous_residual_norm_ > Precision(0) &&
        last_residual_norm_ > previous_residual_norm_) {
        ++residual_increase_count_;
    } else {
        residual_increase_count_ = 0;
    }

    if (previous_residual_norm_ <= Precision(0)) {
        return;
    }

    const Precision relative_change =
        std::abs(last_residual_norm_ - previous_residual_norm_) /
        std::max(previous_residual_norm_, std::numeric_limits<Precision>::epsilon());

    if (relative_change <= stagnation_ratio) {
        ++stagnation_count_;
    } else {
        stagnation_count_ = 0;
    }
}

bool NewtonSolver::should_stop_early_() {
    if (!early_failure_detection) {
        return false;
    }

    if (iterations_ < convergence_check_start) {
        return false;
    }

    const Precision safe_initial = std::max(
        initial_residual_norm_,
        std::numeric_limits<Precision>::epsilon()
    );

    if (last_residual_norm_ > divergence_factor * safe_initial) {
        failed_by_divergence_ = true;
        return true;
    }

    if (maximum_residual_increases > 0 &&
        residual_increase_count_ >= maximum_residual_increases) {
        failed_by_residual_increase_ = true;
        return true;
    }

    if (maximum_stagnation_steps > 0 &&
        stagnation_count_ >= maximum_stagnation_steps) {
        failed_by_stagnation_ = true;
        return true;
    }

    const Precision reduction = last_residual_norm_ / safe_initial;

    if (minimum_residual_reduction > Precision(0) &&
        reduction > minimum_residual_reduction) {
        failed_by_poor_reduction_ = true;
        return true;
    }

    return false;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
