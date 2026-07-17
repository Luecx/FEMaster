/**
 * @file arc_length_control.cpp
 * @brief Implements a generic nonlinear arc-length solver.
 */

#include "arc_length_control.h"

#include "../../core/logging.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {
namespace loadcase {
namespace tools {

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
    logging::error(maximum_increments > 0,
        "ArcLengthControl requires maximum_increments > 0");
    logging::error(maximum_iterations > 0,
        "ArcLengthControl requires maximum_iterations > 0");
    logging::error(tolerance > Precision(0),
        "ArcLengthControl requires tolerance > 0");
    logging::error(initial_increment > Precision(0),
        "ArcLengthControl requires initial_increment > 0");
    logging::error(minimum_increment > Precision(0),
        "ArcLengthControl requires minimum_increment > 0");
    logging::error(maximum_increment >= minimum_increment,
        "ArcLengthControl requires maximum_increment >= minimum_increment");
    logging::error(initial_increment >= minimum_increment &&
                   initial_increment <= maximum_increment,
        "ArcLengthControl requires initial_increment between minimum_increment and maximum_increment");
    logging::error(psi >= Precision(0),
        "ArcLengthControl requires psi >= 0");
    logging::error(growth_factor > Precision(0),
        "ArcLengthControl requires growth_factor > 0");
    logging::error(cutback_factor > Precision(0) && cutback_factor < Precision(1),
        "ArcLengthControl requires cutback_factor between 0 and 1");
    logging::error(fast_iterations > 0,
        "ArcLengthControl requires fast_iterations > 0");
    logging::error(slow_iterations >= fast_iterations,
        "ArcLengthControl requires slow_iterations >= fast_iterations");
    logging::error(maximum_cutbacks > 0,
        "ArcLengthControl requires maximum_cutbacks > 0");

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

    Index cutback_count = 0;

    while (lambda < Precision(1) - tolerance &&
           accepted_increments_ < maximum_increments) {
        q_accepted_      = q;
        lambda_accepted_ = lambda;

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

        if (accepted_increments_ > 0) {
            const Precision direction =
                previous_delta_q_.dot(dq_load)
              + psi2 * load_scale2_ * previous_delta_lambda_;

            path_sign = direction >= Precision(0) ? Precision(1) : Precision(-1);
        }

        arc_radius_ = increment_ * radius_scale_;

        Precision predictor_delta_lambda =
            path_sign * arc_radius_ / predictor_norm;

        const Precision remaining_delta_lambda = Precision(1) - lambda_accepted_;

        if (path_sign > Precision(0) &&
            !adjusting_final_step_ &&
            predictor_delta_lambda >= remaining_delta_lambda) {
            arc_radius_            = remaining_delta_lambda * predictor_norm;
            increment_             = arc_radius_ / radius_scale_;
            predictor_delta_lambda = remaining_delta_lambda;
            adjusting_final_step_  = true;
        }

        if (increment_ < minimum_increment) {
            failure_reason_ = "MINIMUM_INCREMENT";
            return false;
        }

        q      = q_accepted_ + predictor_delta_lambda * dq_load;
        lambda = lambda_accepted_ + predictor_delta_lambda;

        const Precision target_lambda = lambda;

        configure_newton_();

        Precision current_equilibrium_norm = Precision(0);

        bool        converged              = false;
        const char* attempt_failure_reason = "NONE";

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
                logging::error(reference_load.size() == n,
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

                for (Index col = 0; col < tangent.outerSize(); ++col) {
                    for (SparseMatrix::InnerIterator it(tangent, col); it; ++it) {
                        triplets.emplace_back(it.row(), it.col(), it.value());
                    }
                }

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

                logging::error(solution.rows() == n + 1 &&
                               solution.cols() == 1,
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
                Precision convergence_order,
                Time      assembly_ms,
                Time      solve_ms,
                bool      iteration_converged) {
                (void) current_residual_norm;

                if (on_iteration) {
                    on_iteration(
                        accepted_increments_ + 1,
                        iteration,
                        lambda,
                        current_equilibrium_norm,
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
            attempt_failure_reason = newton_.failure_reason();
        }

        if (!converged) {
            q      = q_accepted_;
            lambda = lambda_accepted_;

            if (!adaptive) {
                failure_reason_ = "FIXED_INCREMENT_FAILED";
                return false;
            }

            increment_ *= cutback_factor;
            ++cutback_count;

            logging::info(true,
                "Increment rejected at lambda = ",
                target_lambda,
                "; reason = ",
                attempt_failure_reason,
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
                failure_reason_ = "INVALID_FINAL_DIRECTION";
                return false;
            }

            const Precision previous_increment = increment_;
            arc_radius_ *= remaining_delta_lambda / achieved_delta_lambda;
            increment_   = arc_radius_ / radius_scale_;

            q      = q_accepted_;
            lambda = lambda_accepted_;

            ++cutback_count;

            logging::info(true,
                "Arc-length step reached lambda = ",
                reached_lambda,
                " instead of 1; adjusting increment from ",
                previous_increment,
                " to ",
                increment_
            );

            if (increment_ < minimum_increment) {
                failure_reason_ = "MINIMUM_INCREMENT";
                return false;
            }

            if (cutback_count > maximum_cutbacks) {
                failure_reason_ = "MAXIMUM_FINAL_ADJUSTMENTS";
                return false;
            }

            continue;
        }

        previous_delta_q_      = q - q_accepted_;
        previous_delta_lambda_ = lambda - lambda_accepted_;
        adjusting_final_step_  = false;

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

Index ArcLengthControl::accepted_increments() const {
    return accepted_increments_;
}

Precision ArcLengthControl::increment() const {
    return increment_;
}

const char* ArcLengthControl::failure_reason() const {
    return failure_reason_;
}

void ArcLengthControl::reset_state_() {
    accepted_increments_        = 0;
    increment_                  = initial_increment;
    q_accepted_.resize(0);
    lambda_accepted_            = Precision(0);
    previous_delta_q_.resize(0);
    previous_delta_lambda_      = Precision(1);
    load_scale2_                = Precision(0);
    radius_scale_               = Precision(0);
    arc_radius_                 = Precision(0);
    adjusting_final_step_       = false;
    failure_reason_             = "NONE";
}

void ArcLengthControl::configure_newton_() {
    newton_.maximum_iterations       = maximum_iterations;
    newton_.residual_tolerance       = tolerance;
    newton_.stagnation_tolerance     = Precision(1e-3) * tolerance;
    newton_.check_finite             = true;
    newton_.early_failure_detection  = true;
}

void ArcLengthControl::adapt_increment_() {
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
