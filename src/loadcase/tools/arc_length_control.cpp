/**
 * @file arc_length_control.cpp
 * @brief Implements nonlinear arc-length control.
 *
 * The implementation follows nonlinear static equilibrium paths by solving the
 * finite-element residual together with a scalar arc-length constraint. The path
 * radius is adapted between accepted increments, while each attempted increment
 * may be rolled back and retried when Newton convergence or final-load
 * adjustment fails.
 *
 * All finite-element-specific work is delegated to callbacks supplied by the
 * owning load case. This file only handles the reduced path-following algorithm,
 * augmented Newton system, cutback/growth decisions, active-set restarts and
 * state transaction wiring for nonlinear subsystems.
 *
 * @see ArcLengthControl
 * @see NewtonSolver
 */

#include "arc_length_control.h"

#include "../../core/logging.h"
#include "../../mattools/append_sparse_block.h"

#include <algorithm>
#include <cmath>
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
 * Solves the nonlinear equilibrium path with an arc-length constraint.
 *
 * Each attempted increment starts from the last accepted reduced state `q` and
 * load factor `lambda`. The first tangent predictor determines a load-direction
 * correction and the current radius. Subsequent Newton iterations solve an
 * augmented system containing both the reduced equilibrium residual and the
 * spherical path constraint
 *
 *     ||Delta q||^2 + psi^2 load_scale^2 Delta lambda^2 = radius^2.
 *
 * Increment trial callbacks wrap the complete attempted path increment. They are
 * committed only after the increment is accepted and rolled back for cutbacks,
 * failed Newton solves, failed final-load adjustments or callback exceptions.
 *
 * The current arc-length implementation disables Newton line search because the
 * generic line search scales only the explicit state vector while arc-length
 * corrections also contain an implicit load-factor correction.
 *
 * After a converged augmented Newton solve, `update_active_set` may update
 * discontinuous nonlinear state at the converged configuration. A reported
 * change restarts Newton on the same path increment using the new state.
 *
 * @param q Reduced nonlinear unknown vector. On success it contains the final
 *          accepted state. On failure it is restored to the last accepted state.
 * @param lambda Load factor associated with `q`.
 * @param reference_load Reduced reference load vector used by the arc-length
 *                       predictor and augmented Newton system.
 * @param evaluate Residual and tangent assembly at a supplied `q` and load
 *                 factor.
 * @param linear_solve Linearized solve with one right-hand side.
 * @param matrix_solve Linearized solve for the augmented system.
 * @param residual_norm Problem-specific equilibrium-residual measure.
 * @param correction_norm Problem-specific correction-size measure.
 * @param on_iteration Optional per-Newton-iteration reporting callback.
 * @param on_increment Optional accepted-increment callback for result writing.
 *
 * @return `true` when the path reaches load factor one, otherwise `false`.
 */
bool ArcLengthControl::solve(
    DynamicVector&           q,
    Precision&               lambda,
    const DynamicVector&     reference_load,
    const Evaluate&          evaluate,
    const LinearSolve&       linear_solve,
    const MatrixSolve&       matrix_solve,
    const ResidualNorm&      residual_norm,
    const CorrectionNorm&    correction_norm,
    const IterationCallback& on_iteration,
    const IncrementCallback& on_increment
) {
    logging::error(maximum_iterations > 0,
        "ArcLengthControl requires maximum_iterations > 0");
    logging::error(tolerance > Precision(0),
        "ArcLengthControl requires tolerance > 0");
    logging::error(psi >= Precision(0),
        "ArcLengthControl requires psi >= 0");

    logging::error(static_cast<bool>(evaluate),
        "ArcLengthControl requires an evaluate callback");
    logging::error(static_cast<bool>(linear_solve),
        "ArcLengthControl requires a linear solve callback");
    logging::error(static_cast<bool>(matrix_solve),
        "ArcLengthControl requires a matrix solve callback");
    logging::error(static_cast<bool>(residual_norm),
        "ArcLengthControl requires a residual norm callback");
    logging::error(static_cast<bool>(correction_norm),
        "ArcLengthControl requires a correction norm callback");

    reset_state_();

    previous_delta_q_ = DynamicVector::Zero(q.size());

    while (lambda < Precision(1) - tolerance &&
           increment_manager.can_attempt()) {
        q_accepted_      = q;
        lambda_accepted_ = lambda;

        if (increment_manager.increment() < increment_manager.minimum_increment) {
            increment_manager.fail("MINIMUM_INCREMENT");
            return false;
        }

        Precision current_equilibrium_norm = Precision(0);

        bool        converged              = false;
        const char* attempt_failure_reason = "NONE";
        std::string attempt_error_message;
        Index       active_set_updates     = 0;

        // Treat analysis failures inside this attempted increment like Newton
        // non-convergence so adaptive arc-length control can reduce the radius.
        try {
            if (begin_increment_trial) {
                begin_increment_trial();
            }

            DynamicVector predictor_residual;
            SparseMatrix  predictor_tangent;

            evaluate(
                q_accepted_,
                lambda_accepted_,
                predictor_residual,
                predictor_tangent
            );

            DynamicVector dq_load = linear_solve(
                predictor_tangent,
                reference_load
            );

            logging::error(dq_load.allFinite(),
                "ArcLengthControl: predictor contains NaN/Inf entries");

            const Precision current_load_scale2 = dq_load.squaredNorm();

            logging::error(current_load_scale2 > Precision(0),
                "ArcLengthControl: zero load predictor");

            const Precision psi2 = psi * psi;

            if (radius_scale_ <= Precision(0)) {
                load_scale2_  = current_load_scale2;
                radius_scale_ = std::sqrt(
                    current_load_scale2 + psi2 * load_scale2_
                );
            }

            const Precision predictor_norm = std::sqrt(
                current_load_scale2 + psi2 * load_scale2_
            );

            Precision path_sign = Precision(1);

            if (increment_manager.accepted_increments() > 0) {
                const Precision direction =
                    previous_delta_q_.dot(dq_load)
                  + psi2 * load_scale2_ * previous_delta_lambda_;

                path_sign = direction >= Precision(0) ? Precision(1) : Precision(-1);
            }

            arc_radius_ = increment_manager.increment() * radius_scale_;

            Precision predictor_delta_lambda =
                path_sign * arc_radius_ / predictor_norm;

            const Precision remaining_delta_lambda = Precision(1) - lambda_accepted_;

            if (path_sign > Precision(0) &&
                !adjusting_final_step_ &&
                predictor_delta_lambda >= remaining_delta_lambda) {
                arc_radius_            = remaining_delta_lambda * predictor_norm;
                increment_manager.set_increment(arc_radius_ / radius_scale_);
                predictor_delta_lambda = remaining_delta_lambda;
                adjusting_final_step_  = true;
            }

            q      = q_accepted_ + predictor_delta_lambda * dq_load;
            lambda = lambda_accepted_ + predictor_delta_lambda;

            while (true) {
                configure_newton_();

                converged = newton_.solve(
                    q,
                    [&](const DynamicVector& current_q,
                        DynamicVector&       residual,
                        SparseMatrix&        tangent) {
                        evaluate(current_q, lambda, residual, tangent);
                    },
                    [&](const SparseMatrix&  tangent,
                        const DynamicVector& residual) {
                        const Index n = tangent.rows();

                        logging::error(tangent.rows() == tangent.cols(),
                            "ArcLengthControl: tangent matrix must be square");
                        logging::error(static_cast<Index>(reference_load.size()) == n,
                            "ArcLengthControl: reference load has invalid size");

                        const DynamicVector delta_q      = q - q_accepted_;
                        const Precision     delta_lambda = lambda - lambda_accepted_;
                        const Precision     constraint   =
                            delta_q.dot(delta_q)
                          + psi2 * load_scale2_ * delta_lambda * delta_lambda
                          - arc_radius_ * arc_radius_;

                        SparseMatrix augmented(n + 1, n + 1);
                        TripletList  triplets;
                        triplets.reserve(
                            static_cast<std::size_t>(tangent.nonZeros())
                          + static_cast<std::size_t>(2 * n)
                          + 1
                        );

                        // Embed the equilibrium tangent before adding load and arc-constraint couplings.
                        mattools::append_sparse_block(triplets, tangent);

                        for (Index i = 0; i < n; ++i) {
                            if (reference_load(i) != Precision(0)) {
                                triplets.emplace_back(i, n, -reference_load(i));
                            }

                            if (delta_q(i) != Precision(0)) {
                                triplets.emplace_back(n, i, Precision(2) * delta_q(i));
                            }
                        }

                        triplets.emplace_back(
                            n,
                            n,
                            Precision(2) * psi2 * load_scale2_ * delta_lambda
                        );

                        augmented.setFromTriplets(triplets.begin(), triplets.end());

                        DynamicMatrix rhs(n + 1, 1);
                        rhs.col(0).head(n) = residual;
                        rhs(n, 0)          = -constraint;

                        const DynamicMatrix solution = matrix_solve(augmented, rhs);

                        logging::error(static_cast<Index>(solution.rows()) == n + 1 && solution.cols() == 1,
                            "ArcLengthControl: augmented solve returned an invalid shape");
                        logging::error(solution.allFinite(),
                            "ArcLengthControl: correction contains NaN/Inf entries");

                        const DynamicVector dq      = solution.col(0).head(n);
                        const Precision     dlambda = solution(n, 0);

                        lambda += dlambda;

                        return dq;
                    },
                    [&](const DynamicVector& residual) {
                        current_equilibrium_norm = residual_norm(residual, lambda);

                        return std::max(
                            current_equilibrium_norm,
                            arc_constraint_norm_(q, lambda)
                        );
                    },
                    correction_norm,
                    [&](Index     iteration,
                        Precision current_residual_norm,
                        Precision current_correction_norm,
                        Index     line_search_iterations,
                        Time      assembly_ms,
                        Time      solve_ms) {
                        (void) current_residual_norm;

                        if (on_iteration) {
                            on_iteration(
                                increment_manager.accepted_increments() + 1,
                                iteration,
                                lambda,
                                current_equilibrium_norm,
                                current_correction_norm,
                                line_search_iterations,
                                assembly_ms,
                                solve_ms
                            );
                        }
                    }
                );

                if (!converged) {
                    attempt_failure_reason = newton_.failure_reason();
                    break;
                }

                if (!update_active_set ||
                    update_active_set(q, lambda)) {
                    break;
                }

                ++active_set_updates;

                logging::info(
                    true,
                    "Nonlinear active set changed; restarting Newton at lambda = ",
                    lambda
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
            const Precision rejected_lambda = lambda;

            q      = q_accepted_;
            lambda = lambda_accepted_;

            if (rollback_increment_trial) {
                rollback_increment_trial();
            }

            // Reduce the current path radius without consuming a successful
            // increment. The manager validates fixed/adaptive policy and
            // stores the cause if the requested cutback cannot be scheduled.
            if (!increment_manager.reject(
                    increment_manager.increment(),
                    false,
                    attempt_error_message
                )) {
                return false;
            }

            logging::info(true,
                "Increment rejected at lambda = ",
                rejected_lambda,
                "; reason = ",
                attempt_failure_reason,
                "; reducing increment to ",
                increment_manager.increment()
            );

            continue;
        }

        if (adjusting_final_step_ &&
            std::abs(lambda - Precision(1)) > tolerance) {
            const Precision reached_lambda         = lambda;
            const Precision achieved_delta_lambda =
                lambda - lambda_accepted_;
            const Precision remaining_delta_lambda =
                Precision(1) - lambda_accepted_;

            if (achieved_delta_lambda <= Precision(0)) {
                q               = q_accepted_;
                lambda          = lambda_accepted_;
                increment_manager.fail("INVALID_FINAL_DIRECTION");

                if (rollback_increment_trial) {
                    rollback_increment_trial();
                }

                return false;
            }

            const Precision previous_increment = increment_manager.increment();
            arc_radius_ *= remaining_delta_lambda / achieved_delta_lambda;
            const Precision adjusted_increment = arc_radius_ / radius_scale_;

            q      = q_accepted_;
            lambda = lambda_accepted_;

            if (rollback_increment_trial) {
                rollback_increment_trial();
            }

            // This is a geometric adjustment after Newton convergence, not
            // a Newton failure. It still consumes a retry allowance and
            // restores the accepted model state before another path attempt.
            if (!increment_manager.retry_with(
                    adjusted_increment,
                    "MAXIMUM_FINAL_ADJUSTMENTS"
                )) {
                return false;
            }

            logging::info(true,
                "Arc-length step reached lambda = ",
                reached_lambda,
                " instead of 1; adjusting increment from ",
                previous_increment,
                " to ",
                increment_manager.increment()
            );

            continue;
        }

        previous_delta_q_      = q - q_accepted_;
        previous_delta_lambda_ = lambda - lambda_accepted_;
        adjusting_final_step_  = false;

        if (commit_increment_trial) {
            commit_increment_trial();
        }

        // Store the accepted increment count and adapt only after committing
        // the converged nonlinear trial state.
        increment_manager.accept(newton_.iterations());

        if (on_increment) {
            on_increment(increment_manager.accepted_increments(), q, lambda);
        }

        logging::info(true,
            "Accepted increment "   , increment_manager.accepted_increments(),
            ": lambda = "           , lambda,
            ", Newton iterations = ", newton_.iterations(),
            ", next increment = "  , increment_manager.increment()
        );
    }

    if (lambda < Precision(1) - tolerance) {
        increment_manager.fail("MAXIMUM_INCREMENTS");
        return false;
    }

    return true;
}

const char* ArcLengthControl::failure_reason() const {
    return increment_manager.failure_reason();
}

void ArcLengthControl::reset_state_() {
    increment_manager.reset();
    q_accepted_.resize(0);
    lambda_accepted_            = Precision(0);
    previous_delta_q_.resize(0);
    previous_delta_lambda_      = Precision(1);
    load_scale2_                = Precision(0);
    radius_scale_               = Precision(0);
    arc_radius_                 = Precision(0);
    adjusting_final_step_       = false;
}

void ArcLengthControl::configure_newton_() {
    newton_.maximum_iterations      = maximum_iterations;
    newton_.residual_tolerance      = tolerance;
    newton_.correction_tolerance    = Precision(1e-2);
    newton_.check_finite            = true;
    newton_.early_failure_detection = true;

    // Arc-length Newton corrections contain both dq and dlambda, but the
    // generic Newton line search can scale only the explicit state vector q.
    // Scaling q while lambda already contains the full correction evaluates a
    // point away from the augmented arc-length linearization, so globalization
    // is left to the adaptive arc-length increment control.
    newton_.line_search_enabled = false;
}

Precision ArcLengthControl::arc_constraint_norm_(
    const DynamicVector& q,
    Precision            lambda
) const {
    const DynamicVector delta_q = q - q_accepted_;
    const Precision delta_lambda = lambda - lambda_accepted_;
    const Precision psi2         = psi * psi;

    const Precision constraint =
        delta_q.dot(delta_q)
      + psi2 * load_scale2_ * delta_lambda * delta_lambda
      - arc_radius_ * arc_radius_;

    return std::abs(constraint) / std::max(
        arc_radius_ * arc_radius_,
        std::numeric_limits<Precision>::epsilon()
    );
}

} // namespace tools
} // namespace loadcase
} // namespace fem