/**
 * @file newton_solver.cpp
 * @brief Implements the generic callback-based Newton solver.
 *
 * The implementation performs residual and tangent evaluation, solution of the
 * linearized Newton system, optional backtracking line search, convergence
 * reporting and configurable early-failure detection.
 *
 * All finite-element-specific operations are delegated to callbacks. The
 * Newton solver therefore does not know how residuals and tangents are
 * assembled, how constraints are represented or which sparse factorization is
 * used.
 *
 * The line search evaluates temporary trial states and supports transactional
 * state callbacks for history-dependent models. A trial can be committed when
 * accepted or rolled back when rejected without exposing the corresponding
 * state-management details to the Newton algorithm.
 *
 * @see NewtonSolver
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#include "newton_solver.h"

#include "../../core/logging.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * Solves a nonlinear algebraic system using Newton iterations.
 *
 * Every iteration consists of the following phases:
 *
 * 1. Evaluate the residual vector and tangent matrix at the current state.
 * 2. Compute the caller-defined residual norm.
 * 3. Update failure-detection counters.
 * 4. Accept the current state when both residual and correction criteria pass.
 * 5. Solve the linearized system for the Newton correction.
 * 6. Optionally perform a backtracking line search.
 * 7. Compute the norm of the accepted correction.
 * 8. Report the iteration through the optional callback.
 * 9. Update the nonlinear state.
 *
 * The linear-solve callback must return a correction compatible with the update
 *
 *     x <- x + dx.
 *
 * The solver does not negate the residual or correction internally.
 *
 * When line search is enabled, trial states are evaluated according to
 *
 *     x_trial = x + alpha dx.
 *
 * Starting with `alpha = 1`, the step length is multiplied by
 * `line_search_reduction` after every rejected trial. A trial satisfies the
 * sufficient-decrease criterion when
 *
 *     ||R(x_trial)|| <= (1 - c alpha) ||R(x)||,
 *
 * where `c` is `line_search_sufficient_decrease`.
 *
 * For non-smooth systems, especially active-set contact formulations, this
 * sufficient-decrease condition can reject every trial even though one trial
 * strictly reduces the same residual norm used for convergence. In that case,
 * the best strictly decreasing trial is evaluated again and accepted.
 *
 * Trial evaluations may modify external nonlinear subsystem state. The optional
 * begin, commit and rollback callbacks provide transactional handling of this
 * state. A rejected trial is rolled back before the next step length is tested.
 * An accepted trial is committed before the nonlinear state vector is updated.
 *
 * A return value of `false` indicates one of the classified failure conditions
 * exposed through `failure_reason`. Exceptions thrown by the supplied
 * callbacks are not consumed by this solver, except during individual
 * line-search trials. A failed line-search trial is interpreted as an invalid
 * trial and assigned an infinite residual norm.
 *
 * @param x Nonlinear state vector. The vector is updated in place after every
 *          accepted Newton correction.
 * @param evaluate Callback assembling residual and tangent at a supplied state.
 * @param linear_solve Callback solving the linearized Newton system.
 * @param residual_norm Callback computing the convergence norm of the residual.
 * @param correction_norm Callback computing the norm of an accepted correction.
 * @param on_iteration Optional reporting callback invoked after every completed
 *                     Newton iteration and for convergence detected before a
 *                     new linear solve.
 * @param evaluate_residual Optional residual-only callback used for
 *                          line-search trials. If it is absent, line search
 *                          falls back to the full residual-and-tangent
 *                          evaluation callback.
 *
 * @return `true` when both convergence tolerances are satisfied, otherwise
 *         `false`.
 */
bool NewtonSolver::solve(
    DynamicVector&           x,
    const Evaluate&          evaluate,
    const LinearSolve&       linear_solve,
    const Norm&              residual_norm,
    const CorrectionNorm&    correction_norm,
    const IterationCallback& on_iteration,
    const EvaluateResidual&  evaluate_residual
) {
    // Validate the Newton iteration and convergence settings
    logging::error(maximum_iterations > 0,
        "NewtonSolver requires maximum_iterations > 0");
    logging::error(residual_tolerance > Precision(0),
        "NewtonSolver requires residual_tolerance > 0");
    logging::error(correction_tolerance > Precision(0),
        "NewtonSolver requires correction_tolerance > 0");
    logging::error(convergence_check_start > 0,
        "NewtonSolver requires convergence_check_start > 0");
    logging::error(slow_convergence_check_start >= convergence_check_start,
        "NewtonSolver requires slow_convergence_check_start >= convergence_check_start");
    logging::error(maximum_slow_convergence_checks > 0,
        "NewtonSolver requires maximum_slow_convergence_checks > 0");
    logging::error(divergence_factor > Precision(1),
        "NewtonSolver requires divergence_factor > 1");
    logging::error(residual_growth_tolerance >= Precision(0),
        "NewtonSolver requires residual_growth_tolerance >= 0");
    logging::error(strong_damping_step_length > Precision(0) &&
                   strong_damping_step_length <= Precision(1),
        "NewtonSolver requires strong_damping_step_length between 0 and 1");
    logging::error(strong_damping_residual_ratio > Precision(0) &&
                   strong_damping_residual_ratio < Precision(1),
        "NewtonSolver requires strong_damping_residual_ratio between 0 and 1");
    logging::error(maximum_strong_damping_steps > 0,
        "NewtonSolver requires maximum_strong_damping_steps > 0");

    // Validate the backtracking line-search settings
    logging::error(maximum_line_search_iterations > 0,
        "NewtonSolver requires maximum_line_search_iterations > 0");
    logging::error(line_search_reduction > Precision(0) && line_search_reduction < Precision(1),
        "NewtonSolver requires line_search_reduction between 0 and 1");
    logging::error(line_search_sufficient_decrease >= Precision(0) && line_search_sufficient_decrease < Precision(0.5),
        "NewtonSolver requires line_search_sufficient_decrease between 0 and 0.5");
    logging::error(minimum_step_length > Precision(0) && minimum_step_length <= Precision(1),
        "NewtonSolver requires minimum_step_length between 0 and 1");

    // Validate all mandatory problem-specific callbacks
    logging::error(static_cast<bool>(evaluate),
        "NewtonSolver requires an evaluate callback");
    logging::error(static_cast<bool>(linear_solve),
        "NewtonSolver requires a linear solve callback");
    logging::error(static_cast<bool>(residual_norm),
        "NewtonSolver requires a residual norm callback");
    logging::error(static_cast<bool>(correction_norm),
        "NewtonSolver requires a correction norm callback");

    // Remove all diagnostics and failure flags left by a previous nonlinear
    // solve using the same solver object
    reset_state_();

    // Reuse the residual and tangent objects across Newton iterations to allow
    // the supplied evaluation callback to reuse their allocated storage
    DynamicVector residual{};
    SparseMatrix  tangent{};

    // Perform at most the configured number of Newton iterations
    for (Index iteration = 1; iteration <= maximum_iterations; ++iteration) {
        Timer assembly_timer;
        Timer solve_timer;

        // Assemble the residual and consistent or algorithmic tangent at the
        // current nonlinear state
        assembly_timer.start();
        evaluate(x, residual, tangent);
        assembly_timer.stop();

        iterations_ = iteration;

        // Shift the previously evaluated residual before storing the current one
        update_residual_history_();

        // Evaluate the problem-specific residual convergence measure
        last_residual_norm_ = residual_norm(residual);

        logging::error(!check_finite || std::isfinite(last_residual_norm_),
            "NewtonSolver: residual norm is NaN or Inf");

        // Store the first residual as the reference for hard-divergence checks
        // and initialize the best residual reached in this Newton solve.
        if (iteration == 1) {
            initial_residual_norm_ = last_residual_norm_;
            best_residual_norm_    = last_residual_norm_;
        } else {
            update_failure_state_();
        }

        // Equilibrium and correction convergence are checked before another
        // linear solve. The first evaluation has no preceding correction, so it
        // may converge without the correction criterion only for a practically
        // exact residual. This preserves the linear-response shortcut without
        // accepting a merely force-converged nonlinear state.
        const bool residual_converged =
            last_residual_norm_ <= residual_tolerance;
        const bool first_iteration_converged =
            iteration == 1 &&
            last_residual_norm_ <= std::min(residual_tolerance, Precision(1e-8));
        const bool correction_converged =
            iteration > 1 && last_correction_norm_ <= correction_tolerance;

        if (first_iteration_converged ||
            (residual_converged && correction_converged)) {
            last_step_length_ = Precision(0);

            if (on_iteration) {
                on_iteration(
                    iteration,
                    last_residual_norm_,
                    last_correction_norm_,
                    Index(0),
                    assembly_timer.elapsed(),
                    Time(0)
                );
            }

            return true;
        }

        // Stop before solving another linear system when the residual history
        // already satisfies one of the configured early-failure conditions
        if (should_stop_early_()) {
            if (on_iteration) {
                on_iteration(
                    iteration,
                    last_residual_norm_,
                    last_correction_norm_,
                    Index(0),
                    assembly_timer.elapsed(),
                    Time(0)
                );
            }

            return false;
        }

        // Solve the caller-defined linearized Newton system
        solve_timer.start();
        const DynamicVector full_correction = linear_solve(tangent, residual);
        solve_timer.stop();

        logging::error(!check_finite || full_correction.allFinite(),
            "NewtonSolver: correction contains NaN or Inf entries");

        // Without line search, the complete Newton correction is accepted
        // immediately
        Precision     accepted_step_length   = Precision(1);
        DynamicVector accepted_correction    = full_correction;
        bool          step_accepted          = !line_search_enabled;
        Index         line_search_iterations = 0;

        // Reduce the Newton correction when the full step does not provide
        // sufficient residual decrease
        if (line_search_enabled) {
            Precision trial_step_length = Precision(1);

            // Preserve the best finite trial in case no step satisfies the
            // formal sufficient-decrease condition
            Precision best_step_length = Precision(0);
            Precision best_residual    = std::numeric_limits<Precision>::infinity();

            // Trial acceptance only needs the residual norm. A tangent object
            // is kept for the backwards-compatible path where no residual-only
            // callback has been supplied by the nonlinear strategy.
            DynamicVector trial_residual{};
            SparseMatrix  trial_tangent{};

            // Indicates whether an external transactional trial state has been
            // opened and still requires commit or rollback
            bool trial_open = false;

            // Begin a temporary external state for one trial evaluation
            const auto begin_trial = [&]() {
                if (begin_line_search_trial) {
                    begin_line_search_trial();
                }
            };

            // Accept the external state produced by the current trial
            const auto commit_trial = [&]() {
                if (commit_line_search_trial) {
                    commit_line_search_trial();
                }
            };

            // Restore the external state that existed before the current trial
            const auto rollback_trial = [&]() {
                if (rollback_line_search_trial) {
                    rollback_line_search_trial();
                }
            };

            /**
             * Evaluates one line-search trial and returns its residual norm.
             *
             * Any exception raised during trial-state preparation, residual
             * assembly or norm evaluation marks the trial as invalid. If the
             * temporary state was opened successfully, it is rolled back before
             * the infinite trial norm is returned.
             */
            const auto evaluate_trial = [&](Precision step_length) {
                const DynamicVector trial_state = x + step_length * full_correction;

                try {
                    // Open a fresh trial transaction before evaluating the
                    // temporary nonlinear state
                    begin_trial();
                    trial_open = true;

                    if (evaluate_residual) {
                        evaluate_residual(trial_state, trial_residual);
                    } else {
                        evaluate(trial_state, trial_residual, trial_tangent);
                    }

                    // Invalid residual entries are treated like a failed trial
                    // without passing them to an arbitrary norm implementation
                    if (!trial_residual.allFinite()) {
                        return std::numeric_limits<Precision>::infinity();
                    }

                    const Precision trial_norm = residual_norm(trial_residual);

                    return std::isfinite(trial_norm)
                        ? trial_norm
                        : std::numeric_limits<Precision>::infinity();
                } catch (...) {
                    // Restore external nonlinear subsystem state before
                    // classifying the trial as invalid
                    if (trial_open) {
                        rollback_trial();
                        trial_open = false;
                    }

                    return std::numeric_limits<Precision>::infinity();
                }
            };

            // Test successively smaller corrections until a trial is accepted
            // or the line-search limits are reached
            for (Index line_search_iteration = 0;
                 line_search_iteration < maximum_line_search_iterations;
                 ++line_search_iteration) {
                const Precision trial_norm = evaluate_trial(trial_step_length);
                ++line_search_iterations;

                // Remember the best finite trial independently of the formal
                // sufficient-decrease condition
                if (std::isfinite(trial_norm) && trial_norm < best_residual) {
                    best_residual    = trial_norm;
                    best_step_length = trial_step_length;
                }

                // Armijo-like sufficient decrease based on the same residual
                // norm used for nonlinear convergence
                const Precision required_norm =
                    (Precision(1) - line_search_sufficient_decrease * trial_step_length) *
                    last_residual_norm_;

                const bool sufficient_decrease =
                    std::isfinite(trial_norm) &&
                    trial_norm <= required_norm;

                // Preserve the accepted external trial state and the
                // corresponding scaled correction
                if (sufficient_decrease) {
                    accepted_step_length = trial_step_length;
                    accepted_correction  = trial_step_length * full_correction;
                    step_accepted        = true;

                    commit_trial();
                    trial_open = false;

                    break;
                }

                // Every rejected trial must restore the external state before
                // another step length is evaluated
                if (trial_open) {
                    rollback_trial();
                    trial_open = false;
                }

                // Reduce the trial step geometrically
                trial_step_length *= line_search_reduction;

                if (trial_step_length < minimum_step_length) {
                    break;
                }
            }

            // Non-smooth active-set changes can violate the smooth
            // sufficient-decrease model even when one trial genuinely improves
            // the residual. Re-evaluate and accept the best strictly decreasing
            // trial in that case.
            if (!step_accepted &&
                best_step_length > Precision(0) &&
                best_residual < last_residual_norm_) {
                const Precision retry_norm = evaluate_trial(best_step_length);
                ++line_search_iterations;

                if (std::isfinite(retry_norm) && retry_norm < last_residual_norm_) {
                    accepted_step_length = best_step_length;
                    accepted_correction  = best_step_length * full_correction;
                    step_accepted        = true;

                    commit_trial();
                    trial_open = false;
                } else if (trial_open) {
                    rollback_trial();
                    trial_open = false;
                }
            }

            // No admissible line-search step was found. Keep the current state
            // vector unchanged and classify the nonlinear solve as failed.
            if (!step_accepted) {
                last_step_length_      = best_step_length;
                failed_by_line_search_ = true;

                return false;
            }
        }

        // Store the accepted step length and evaluate the problem-specific
        // correction measure before updating the nonlinear state
        last_step_length_     = accepted_step_length;
        last_correction_norm_ = correction_norm(x, accepted_correction);

        logging::error(!check_finite || std::isfinite(last_correction_norm_),
            "NewtonSolver: correction norm is NaN or Inf");

        // Report the completed non-converged iteration before applying the
        // accepted correction
        if (on_iteration) {
            on_iteration(
                iteration,
                last_residual_norm_,
                last_correction_norm_,
                line_search_iterations,
                assembly_timer.elapsed(),
                solve_timer.elapsed()
            );
        }

        // Advance the nonlinear state using the complete or line-search-scaled
        // Newton correction
        x += accepted_correction;

    }

    // The convergence tolerances were not reached within the configured
    // iteration limit
    failed_by_maximum_iterations_ = true;

    return false;
}

/**
 * Returns the number of residual evaluations performed by the most recent
 * nonlinear solve.
 *
 * @return Number of completed Newton iterations.
 */
Index NewtonSolver::iterations() const {
    return iterations_;
}

/**
 * Returns the step length of the most recently accepted Newton correction.
 *
 * A value below one indicates that the line search reduced the full Newton
 * correction.
 *
 * @return Last accepted line-search step length.
 */
Precision NewtonSolver::last_step_length() const {
    return last_step_length_;
}

/**
 * Returns a stable textual classification of the most recent solver failure.
 *
 * The failure conditions are checked in the same priority order in which they
 * are reported. A successful solve or a solver that has not yet run returns
 * `"NONE"`.
 *
 * @return Null-terminated failure-reason identifier.
 */
const char* NewtonSolver::failure_reason() const {
    if (failed_by_divergence_) {
        return "DIVERGENCE";
    }

    if (failed_by_strong_damping_) {
        return "LINE_SEARCH_STAGNATION";
    }

    if (failed_by_slow_convergence_) {
        return "SLOW_CONVERGENCE";
    }

    if (failed_by_line_search_) {
        return "LINE_SEARCH_FAILED";
    }

    if (failed_by_maximum_iterations_) {
        return "MAXIMUM_ITERATIONS";
    }

    return "NONE";
}

/**
 * Resets all runtime diagnostics and failure classification.
 *
 * Configuration values and callback objects remain unchanged so the solver can
 * be reused with the same numerical settings.
 */
void NewtonSolver::reset_state_() {
    // Reset iteration and convergence diagnostics
    iterations_ = 0;

    initial_residual_norm_         = Precision(0);
    best_residual_norm_            = Precision(0);
    last_residual_norm_            = Precision(0);
    last_correction_norm_          = Precision(0);
    previous_residual_norm_        = Precision(0);
    last_step_length_              = Precision(1);
    projected_iteration_count_     = Precision(0);
    strong_damping_reference_norm_ = Precision(0);

    // Reset persistent failure evidence
    slow_convergence_count_ = 0;
    strong_damping_count_   = 0;

    // Reset all mutually descriptive failure flags
    failed_by_divergence_         = false;
    failed_by_slow_convergence_   = false;
    failed_by_strong_damping_     = false;
    failed_by_line_search_        = false;
    failed_by_maximum_iterations_ = false;
}

/**
 * Stores the previously evaluated residual norm before the current norm is
 * assigned.
 */
void NewtonSolver::update_residual_history_() {
    previous_residual_norm_ = last_residual_norm_;
}

/**
 * Updates persistent evidence used by the early-failure detector.
 *
 * The slow-convergence forecast follows the idea used by CalculiX: the ratio
 *
 *     rho = r_i / r_(i-1)
 *
 * is interpreted as the current asymptotic residual decay rate. For 0 < rho < 1,
 * the estimated total iteration count required to reach the residual tolerance is
 *
 *     i_est = i + ceil(log(tol / r_i) / log(rho)).
 *
 * A single pessimistic estimate is not enough to reject the increment because
 * Newton convergence can accelerate abruptly after yielding or active-set changes.
 * Several consecutive forecasts beyond the configured iteration budget are
 * required.
 *
 * Strong line-search damping is tracked independently. Several accepted steps at
 * or below strong_damping_step_length only count as stagnation when the residual
 * has failed to decrease by at least the configured amount across the whole
 * sequence.
 */
void NewtonSolver::update_failure_state_() {
    const Precision epsilon = std::numeric_limits<Precision>::epsilon();

    // CalculiX-inspired convergence-rate forecast.
    if (iterations_ >= slow_convergence_check_start &&
        previous_residual_norm_ > residual_tolerance &&
        last_residual_norm_ > residual_tolerance) {
        const Precision rate =
            last_residual_norm_ /
            std::max(previous_residual_norm_, epsilon);

        bool projected_beyond_budget = false;
        projected_iteration_count_   = Precision(0);

        if (rate >= Precision(1)) {
            projected_iteration_count_ =
                std::numeric_limits<Precision>::infinity();
            projected_beyond_budget = true;
        } else if (rate > Precision(0)) {
            const Precision remaining =
                std::log(residual_tolerance / last_residual_norm_) /
                std::log(rate);

            if (std::isfinite(remaining)) {
                projected_iteration_count_ =
                    Precision(iterations_) + std::ceil(std::max(Precision(0), remaining));
                projected_beyond_budget =
                    projected_iteration_count_ > Precision(maximum_iterations);
            }
        }

        if (projected_beyond_budget) {
            ++slow_convergence_count_;
        } else {
            slow_convergence_count_ = 0;
        }
    } else {
        slow_convergence_count_ = 0;
        projected_iteration_count_ = Precision(0);
    }

    // A strongly damped accepted step is meaningful evidence only when it repeats.
    // Keep the residual at the start of the damping sequence so the group can be
    // judged by actual progress rather than by alpha alone.
    if (last_step_length_ <= strong_damping_step_length) {
        if (strong_damping_count_ == 0) {
            strong_damping_reference_norm_ = previous_residual_norm_;
        }
        ++strong_damping_count_;
    } else {
        strong_damping_count_          = 0;
        strong_damping_reference_norm_ = Precision(0);
    }

    best_residual_norm_ = std::min(best_residual_norm_, last_residual_norm_);
}

/**
 * Evaluates robust early-failure criteria.
 *
 * The checks intentionally distinguish three situations:
 *
 * 1. clear divergence: catastrophic growth or a sustained move away from the
 *    best residual reached so far;
 * 2. line-search stagnation: several severely damped accepted steps that achieve
 *    less than the required aggregate residual reduction;
 * 3. slow convergence: several consecutive CalculiX-style convergence forecasts
 *    predict that the residual tolerance cannot be reached within the Newton
 *    iteration budget.
 *
 * This avoids the former sensitivity to one arbitrary residual-ratio threshold
 * or one tiny correction while still cutting increments that are demonstrably
 * wasting nonlinear iterations.
 */
bool NewtonSolver::should_stop_early_() {
    if (!early_failure_detection) {
        return false;
    }

    const Precision epsilon = std::numeric_limits<Precision>::epsilon();
    const Precision safe_initial =
        std::max(initial_residual_norm_, epsilon);

    // Catastrophic divergence can be detected immediately after the first step.
    if (iterations_ > 1 &&
        last_residual_norm_ > divergence_factor * safe_initial) {
        failed_by_divergence_ = true;
        return true;
    }

    // CalculiX similarly compares the current residual with the best residual
    // already reached before declaring divergence. Add a small relative margin
    // so harmless numerical oscillation around the minimum is not classified as
    // divergence.
    if (iterations_ >= convergence_check_start &&
        best_residual_norm_ > Precision(0)) {
        const Precision growth_limit =
            (Precision(1) + residual_growth_tolerance) * best_residual_norm_;

        if (previous_residual_norm_ > growth_limit &&
            last_residual_norm_ > growth_limit &&
            last_residual_norm_ > previous_residual_norm_ &&
            last_residual_norm_ > residual_tolerance) {
            failed_by_divergence_ = true;
            return true;
        }
    }

    // Repeated severe damping is a failure only when the complete sequence made
    // less than (1 - strong_damping_residual_ratio) relative residual progress.
    if (strong_damping_count_ >= maximum_strong_damping_steps &&
        strong_damping_reference_norm_ > Precision(0) &&
        last_residual_norm_ >
            strong_damping_residual_ratio * strong_damping_reference_norm_) {
        failed_by_strong_damping_ = true;
        return true;
    }

    // Do not trust one convergence forecast: plastic transitions can suddenly
    // accelerate after several mediocre iterations. Persistent forecasts beyond
    // the actual Newton budget are strong evidence that a smaller increment is
    // cheaper than continuing.
    if (slow_convergence_count_ >= maximum_slow_convergence_checks) {
        failed_by_slow_convergence_ = true;
        return true;
    }

    return false;
}

} // namespace tools
} // namespace loadcase
} // namespace fem