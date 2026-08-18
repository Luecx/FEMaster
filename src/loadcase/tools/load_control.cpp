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
 * @param lambda Load factor. On success it reaches one within tolerance. On
 *               failure it is restored to the last accepted value.
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
                                accepted_increments_ + 1,
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

            if (!adaptive) {
                failure_reason_ = attempt_error_message.empty()
                    ? "FIXED_INCREMENT_FAILED"
                    : "FIXED_INCREMENT_FAILED: " + attempt_error_message;
                return false;
            }

            increment_ *= cutback_factor;
            ++cutback_count;

            logging::info(
                true,
                "Increment rejected at lambda = ",
                target_lambda,
                "; reducing increment to ",
                increment_
            );

            if (increment_ < minimum_increment) {
                failure_reason_ = attempt_error_message.empty()
                    ? "MINIMUM_INCREMENT"
                    : "MINIMUM_INCREMENT: " + attempt_error_message;
                return false;
            }

            if (cutback_count > maximum_cutbacks) {
                failure_reason_ = attempt_error_message.empty()
                    ? "MAXIMUM_CUTBACKS"
                    : "MAXIMUM_CUTBACKS: " + attempt_error_message;
                return false;
            }

            continue;
        }

        lambda = target_lambda;

        if (commit_increment_trial) {
            commit_increment_trial();
        }

        ++accepted_increments_;
        cutback_count = 0;

        if (on_increment) {
            on_increment(
                accepted_increments_,
                q,
                lambda
            );
        }

        adapt_increment_();

        logging::info(
            true,
            "Accepted increment ",
            accepted_increments_,
            ": lambda = ",
            lambda,
            ", Newton iterations = ",
            newton_.iterations(),
            ", next increment = ",
            increment_
        );
    }

    if (lambda < Precision(1) - tolerance) {
        failure_reason_ = "MAXIMUM_INCREMENTS";
        return false;
    }

    return true;
}

const char* LoadControl::failure_reason() const {
    return failure_reason_.c_str();
}

void LoadControl::reset_state_() {
    accepted_increments_ = 0;
    increment_           = initial_increment;
    failure_reason_      = "NONE";
}

void LoadControl::configure_newton_() {
    newton_.maximum_iterations      = maximum_iterations;
    newton_.residual_tolerance      = tolerance;
    newton_.stagnation_tolerance    = Precision(1e-3) * tolerance;
    newton_.check_finite            = true;
    newton_.early_failure_detection = false;

    newton_.begin_line_search_trial    = begin_line_search_trial;
    newton_.commit_line_search_trial   = commit_line_search_trial;
    newton_.rollback_line_search_trial = rollback_line_search_trial;
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
