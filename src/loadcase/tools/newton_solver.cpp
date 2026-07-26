/**
 * @file newton_solver.cpp
 * @brief Implements a generic Newton solver.
 */

#include "newton_solver.h"

#include "../../core/logging.h"

#include <algorithm>
#include <cmath>
#include <iostream>
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

    logging::error(maximum_line_search_iterations > 0,
        "NewtonSolver requires maximum_line_search_iterations > 0");
    logging::error(line_search_reduction > Precision(0) &&
                   line_search_reduction < Precision(1),
        "NewtonSolver requires line_search_reduction between 0 and 1");
    logging::error(line_search_sufficient_decrease >= Precision(0) &&
                   line_search_sufficient_decrease < Precision(0.5),
        "NewtonSolver requires line_search_sufficient_decrease between 0 and 0.5");
    logging::error(minimum_step_length > Precision(0) &&
                   minimum_step_length <= Precision(1),
        "NewtonSolver requires minimum_step_length between 0 and 1");

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
        const DynamicVector dx = linear_solve(tangent, residual);
        solve_timer.stop();

        logging::error(
            !check_finite || dx.allFinite(),
            "NewtonSolver: correction contains NaN/Inf entries"
        );

        Precision     accepted_alpha = Precision(1);
        DynamicVector accepted_dx    = dx;
        bool          step_accepted  = !line_search_enabled;

        if (line_search_enabled) {
            Precision alpha = Precision(1);

            Precision best_alpha = Precision(0);
            Precision best_norm  = std::numeric_limits<Precision>::infinity();

            DynamicVector trial_residual;
            SparseMatrix  trial_tangent;

            auto begin_trial = [&]() {
                if (begin_line_search_trial) {
                    begin_line_search_trial();
                }
            };

            auto commit_trial = [&]() {
                if (commit_line_search_trial) {
                    commit_line_search_trial();
                }
            };

            auto rollback_trial = [&]() {
                if (rollback_line_search_trial) {
                    rollback_line_search_trial();
                }
            };

            auto evaluate_trial = [&](Precision trial_alpha) {
                const DynamicVector x_trial = x + trial_alpha * dx;

                begin_trial();

                try {
                    evaluate(x_trial, trial_residual, trial_tangent);

                    Precision trial_norm =
                        std::numeric_limits<Precision>::infinity();

                    if (trial_residual.allFinite()) {
                        trial_norm = residual_norm(trial_residual);
                    }

                    return trial_norm;
                } catch (...) {
                    rollback_trial();
                    throw;
                }
            };

            for (Index line_search_iteration = 0;
                 line_search_iteration < maximum_line_search_iterations;
                 ++line_search_iteration) {
                const Precision trial_norm =
                    evaluate_trial(alpha);

                if (std::isfinite(trial_norm) && trial_norm < best_norm) {
                    best_norm  = trial_norm;
                    best_alpha = alpha;
                }

                const Precision required_norm =
                    (Precision(1) -
                     line_search_sufficient_decrease * alpha) *
                    last_residual_norm_;

                const bool accepted =
                    std::isfinite(trial_norm) &&
                    trial_norm <= required_norm;

                std::cout
                    << "[LINESEARCH]"
                    << " iter="       << iter
                    << " ls="         << line_search_iteration
                    << " alpha="      << alpha
                    << " rel_res="    << last_residual_norm_
                    << " trial_rel="  << trial_norm
                    << " limit="      << required_norm
                    << " accepted="   << accepted
                    << '\n';

                if (accepted) {
                    accepted_alpha = alpha;
                    accepted_dx    = alpha * dx;
                    step_accepted  = true;
                    commit_trial();
                    break;
                }

                rollback_trial();

                alpha *= line_search_reduction;

                if (alpha < minimum_step_length) {
                    break;
                }
            }

            /*
             * For non-smooth active-set changes, the sufficient-decrease
             * condition may be too strict even though a trial genuinely lowers
             * the same residual norm used for convergence. Accept the best
             * strictly decreasing trial in that case.
             */
            if (!step_accepted &&
                best_alpha > Precision(0) &&
                best_norm < last_residual_norm_) {
                (void) evaluate_trial(best_alpha);
                commit_trial();

                accepted_alpha = best_alpha;
                accepted_dx    = best_alpha * dx;
                step_accepted  = true;
            }

            if (!step_accepted) {
                last_step_length_      = best_alpha;
                failed_by_line_search_ = true;

                std::cout
                    << "[LINESEARCH] FAILED"
                    << " iter="      << iter
                    << " best_alpha=" << best_alpha
                    << " best_rel="   << best_norm
                    << '\n';

                return false;
            }
        }

        last_step_length_     = accepted_alpha;
        last_correction_norm_ = correction_norm(x, accepted_dx);

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

        x += accepted_dx;

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

Precision NewtonSolver::last_step_length() const {
    return last_step_length_;
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

bool NewtonSolver::failed_by_line_search() const {
    return failed_by_line_search_;
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

    if (failed_by_line_search_) {
        return "LINE_SEARCH_FAILED";
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
    last_step_length_                = Precision(1);

    residual_increase_count_ = 0;
    stagnation_count_        = 0;

    failed_by_divergence_         = false;
    failed_by_residual_increase_  = false;
    failed_by_stagnation_         = false;
    failed_by_poor_reduction_     = false;
    failed_by_line_search_        = false;
    failed_by_maximum_iterations_ = false;
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

    if (r0 <= Precision(0) ||
        r1 <= Precision(0) ||
        r2 <= Precision(0)) {
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
        std::max(
            previous_residual_norm_,
            std::numeric_limits<Precision>::epsilon()
        );

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

    const Precision reduction =
        last_residual_norm_ / safe_initial;

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
