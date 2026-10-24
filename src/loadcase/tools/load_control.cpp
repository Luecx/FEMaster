/**
 * @file load_control.cpp
 * @brief Implements direct load-factor control for nonlinear equilibrium paths.
 *
 * The implementation advances the scalar load factor by attempted increments,
 * solves the nonlinear residual at each target load factor, and applies adaptive
 * growth or cutback based on Newton convergence. It does not assemble finite
 * element matrices itself; all model-specific operations are supplied by the
 * owning load case through callbacks.
 *
 * Stateful nonlinear subsystems are handled through explicit trial
 * transactions. Increment trials protect state across cutbacks, while
 * line-search trials protect temporary residual evaluations inside the Newton
 * solver. The same mechanism is used by contact and can be reused by nonlinear
 * material history updates.
 *
 * @see LoadControl
 * @see NewtonSolver
 */

#include "load_control.h"

#include "../../core/logging.h"

#include <algorithm>
#include <exception>
#include <limits>
#include <string>

namespace fem {
namespace loadcase {
namespace tools {
namespace {

constexpr Index maximum_active_set_updates = 8;

} // namespace

/**
 * Solves the nonlinear equilibrium path by prescribing the load factor.
 *
 * Each attempted increment starts from the last accepted reduced state `q` and
 * load factor `lambda`. The target load factor is fixed during the Newton solve,
 * so the Newton state contains only the reduced unknowns. A caller-supplied
 * predictor may initialize `q` for the target load factor before the first
 * Newton evaluation.
 *
 * The increment trial callbacks wrap the entire attempted step. They are opened
 * before predictor and Newton evaluations, committed only after the increment is
 * accepted, and rolled back when Newton fails, an analysis callback throws or an
 * adaptive cutback is required.
 *
 * Line-search trial callbacks are forwarded to the Newton solver and are nested
 * inside the active increment trial. They allow temporary residual evaluations
 * to modify trial-only state without contaminating the accepted increment state.
 *
 * After a converged Newton solve, `update_active_set` may update discontinuous
 * nonlinear state at the converged configuration. If it reports a change, Newton
 * is restarted at the same target load factor.
 *
 * @param q Reduced nonlinear unknown vector. On success it contains the final
 *          accepted state. On failure it is restored to the last accepted state.
 * @param lambda Load factor. On success the final equilibrium is solved at
 *               exactly one; on failure the last accepted factor is retained.
 * @param evaluate Residual and tangent assembly at a supplied `q` and load
 *                 factor.
 * @param linear_solve Linearized Newton solve returning the correction added to
 *                     `q`.
 * @param residual_norm Problem-specific residual convergence measure.
 * @param correction_norm Problem-specific correction-size measure.
 * @param on_iteration Optional per-Newton-iteration reporting callback.
 * @param on_increment Optional accepted-increment callback for result writing.
 * @param predictor Optional tangent predictor for the attempted target load.
 * @param evaluate_residual Optional residual-only evaluation used by line
 *                          search.
 *
 * @return `true` when the load factor reaches one, otherwise `false`.
 */
bool LoadControl::solve(
    DynamicVector&           q,
    Precision&               lambda,
    const Evaluate&          evaluate,
    const LinearSolve&       linear_solve,
    const ResidualNorm&      residual_norm,
    const CorrectionNorm&    correction_norm,
    const IterationCallback& on_iteration,
    const IncrementCallback& on_increment,
    const Predictor&         predictor,
    const EvaluateResidual&  evaluate_residual
) {
    logging::error(maximum_iterations > 0,
        "LoadControl requires maximum_iterations > 0");
    logging::error(tolerance > Precision(0),
        "LoadControl requires tolerance > 0");
    logging::error(static_cast<bool>(evaluate),
        "LoadControl requires an evaluate callback");
    logging::error(static_cast<bool>(linear_solve),
        "LoadControl requires a linear solve callback");
    logging::error(static_cast<bool>(residual_norm),
        "LoadControl requires a residual norm callback");
    logging::error(static_cast<bool>(correction_norm),
        "LoadControl requires a correction norm callback");

    reset_state_();
    Precision lambda_roundoff = Precision(0);

    // Newton's residual tolerance is not a permissible error in the final load
    // factor. A step accepted slightly below lambda=1 could leave a nonlinear
    // amplitude at a different value from A(T_step)=1; the next analysis would
    // then inherit the full nominal load and introduce an artificial jump.
    // Always converge an actual equilibrium at lambda=1 before reporting success.
    while (lambda < Precision(1) && increment_manager.accepted_increments() < increment_manager.can_attempt()) {
        const Precision remaining = Precision(1) - lambda;
        const Precision step      = std::min(increment_manager.increment(), remaining);

        // If the proposed endpoint is already within the Newton tolerance,
        // solve at exactly lambda=1 instead of accepting a near-end frame.
        // Repeated additions can leave lambda slightly below its exact value
        // (e.g. 0.8999999999999999 after nine steps of 0.1). Allow only a few
        // machine epsilons when checking the final maximum increment, so this
        // rounding does not prevent the last solve without permitting a
        // genuinely larger load increment. Never snap during cutback retries.
        const Precision endpoint_roundoff = Precision(4) * std::numeric_limits<Precision>::epsilon();
        const bool snap_to_end = increment_manager.cutback_count() == 0
                              && remaining - increment_manager.maximum_increment <= endpoint_roundoff
                              && remaining - step <= tolerance;

        // Compensate for rounding in successive additions of the load factor.
        // Without this correction, many otherwise valid fixed increments can
        // accumulate enough error to miss lambda=1 and exhaust MAX_INCREMENTS.
        // The compensation belongs to accepted increments only: failed trials
        // and cutback retries must start from exactly the same accepted state.
        const Precision corrected_step = step - lambda_roundoff;
        const Precision proposed_lambda = lambda + corrected_step;
        const Precision target_lambda = snap_to_end
            ? Precision(1) : std::min(proposed_lambda, Precision(1));

        // A snap can enlarge the proposed increment. Ordinary cutbacks must
        // retain the original step to avoid cancellation in target-lambda
        // subtraction; snapped cutbacks use the full attempted remainder.
        const Precision attempted_increment = snap_to_end ? remaining : step;

        // Compare the proposed step directly with the configured minimum.
        // Subtracting lambda from target_lambda can lose precision and reject
        // valid fixed increments (e.g. 0.1) due to floating-point roundoff.
        // Only an actual final step to lambda=1 may be below the minimum.
        if (step < increment_manager.minimum_increment
            && target_lambda != Precision(1)) {
            failure_reason_ = "MINIMUM_INCREMENT";
            return false;
        }

        const DynamicVector q_accepted      = q;
        const Precision     lambda_accepted = lambda;

        bool        converged              = false;
        const char* attempt_failure_reason = "NONE";
        std::string attempt_error_message;
        Index       active_set_updates     = 0;

        // Treat analysis failures inside this attempted increment like Newton
        // non-convergence so adaptive load control can reduce the step size.
        try {
            if (begin_increment_trial) {
                begin_increment_trial();
            }

            if (predictor) {
                predictor(q, lambda_accepted, target_lambda);
            }

            while (true) {
                configure_newton_();

                NewtonSolver::EvaluateResidual newton_evaluate_residual;
                if (evaluate_residual) {
                    newton_evaluate_residual =
                        [&](const DynamicVector& current_q,
                            DynamicVector&       residual) {
                            evaluate_residual(
                                current_q,
                                target_lambda,
                                residual
                            );
                        };
                }

                converged = newton_.solve(
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
                        Index     line_search_iterations,
                        Time      assembly_ms,
                        Time      solve_ms) {
                        if (on_iteration) {
                            on_iteration(
                                increment_manager.accepted_increments() + 1,
                                iteration,
                                target_lambda,
                                current_residual_norm,
                                current_correction_norm,
                                line_search_iterations,
                                assembly_ms,
                                solve_ms
                            );
                        }
                    },
                    newton_evaluate_residual
                );

                if (!converged) {
                    attempt_failure_reason = newton_.failure_reason();
                    break;
                }

                if (!update_active_set ||
                    update_active_set(q, target_lambda)) {
                    break;
                }

                ++active_set_updates;

                logging::info(
                    true,
                    "Nonlinear active set changed; restarting Newton at lambda = ",
                    target_lambda
                );

                if (active_set_updates > maximum_active_set_updates) {
                    converged              = false;
                    attempt_failure_reason = "ACTIVE_SET";
                    break;
                }
            }
        } catch (const std::exception& exception) {
            converged              = false;
            attempt_failure_reason = "ANALYSIS_ERROR";
            attempt_error_message  = exception.what();
        } catch (...) {
            converged              = false;
            attempt_failure_reason = "ANALYSIS_ERROR";
            attempt_error_message  = "unknown exception";
        }

        if (!converged) {
            q      = q_accepted;
            lambda = lambda_accepted;

            if (rollback_increment_trial) {
                rollback_increment_trial();
            }

            logging::info(
                true,
                "Newton failed: ",
                attempt_failure_reason,
                ", last alpha = ",
                newton_.last_step_length()
            );

            const std::string failure(attempt_failure_reason);
            const bool severe_failure =
                failure == "DIVERGENCE" ||
                failure == "LINE_SEARCH_FAILED" ||
                failure == "LINE_SEARCH_STAGNATION" ||
                failure == "ANALYSIS_ERROR";

            // The manager applies a stronger cutback for serious Newton
            // failures, checks fixed-increment policy and retry limits, and
            // retains diagnostics without consuming an accepted increment.
            if (!increment_manager.reject(
                    attempted_increment,
                    severe_failure,
                    attempt_error_message
                )) {
                return false;
            }

            logging::info(
                true,
                "Increment rejected at lambda = ",
                target_lambda,
                "; reducing increment to ",
                increment_manager.increment()
            );

            continue;
        }

        // Kahan compensation stores the rounding error in the accepted
        // addition. Do not update it on failed trials, and clear it once the
        // exact endpoint has been reached.
        lambda_roundoff = target_lambda == Precision(1)
            ? Precision(0)
            : (target_lambda - lambda_accepted) - corrected_step;
        lambda = target_lambda;

        if (commit_increment_trial) {
            commit_increment_trial();
        }

        // Record this converged increment and adapt the next nominal step.
        // Rejected trials do not advance the accepted-increment counter.
        increment_manager.accept(newton_.iterations());

        if (on_increment) {
            on_increment(
                increment_manager.accepted_increments(),
                q,
                lambda
            );
        }

        logging::info(
            true,
            "Accepted increment ",
            increment_manager.accepted_increments(),
            ": lambda = ",
            lambda,
            ", Newton iterations = ",
            newton_.iterations(),
            ", next increment = ",
            increment_manager.increment()
        );
    }

    if (lambda < Precision(1)) {
        increment_manager.fail("MAXIMUM_INCREMENTS");
        return false;
    }

    return true;
}

const char* LoadControl::failure_reason() const {
    return increment_manager.failure_reason();
}

void LoadControl::reset_state_() {
    increment_manager.reset();
}

void LoadControl::configure_newton_() {
    newton_.maximum_iterations      = maximum_iterations;
    newton_.residual_tolerance      = tolerance;
    newton_.correction_tolerance    = Precision(1e-2);
    newton_.check_finite            = true;
    newton_.early_failure_detection = adaptive;

    newton_.begin_line_search_trial    = begin_line_search_trial;
    newton_.commit_line_search_trial   = commit_line_search_trial;
    newton_.rollback_line_search_trial = rollback_line_search_trial;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
