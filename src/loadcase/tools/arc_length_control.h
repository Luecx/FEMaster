/**
 * @file arc_length_control.h
 * @brief Declares arc-length control for nonlinear equilibrium path following.
 *
 * Arc-length control advances a nonlinear static analysis by constraining the
 * combined length of the displacement correction and load-factor correction.
 * This allows the nonlinear solver to follow equilibrium paths whose load
 * factor is not monotonic, such as snap-through or post-buckling responses.
 *
 * The class owns only the scalar path constraint, increment radius management,
 * cutback/growth decisions, active-set restarts and callback wiring. Residual
 * assembly, tangent assembly, constraint projection, sparse solves, convergence
 * measures and result writing remain responsibilities of the owning load case.
 *
 * Increment trial callbacks protect stateful nonlinear subsystems across failed
 * path increments and cutbacks.
 *
 * @see ArcLengthControl
 * @see NewtonSolver
 */

#pragma once

#include "newton_solver.h"

#include <functional>
#include <string>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * @brief Nonlinear path controller using a spherical arc-length constraint.
 *
 * `ArcLengthControl` solves the reduced nonlinear equilibrium system together
 * with an additional scalar constraint
 *
 *     ||Delta q||^2 + psi^2 load_scale^2 Delta lambda^2 = radius^2.
 *
 * The radius is initialized from the first load predictor and then adapted
 * between accepted increments. During an increment, Newton solves the augmented
 * system for both the reduced displacement correction and the load-factor
 * correction. The current implementation disables Newton line search because a
 * generic step scaling of `q` alone would be inconsistent with the coupled
 * `Delta lambda` correction.
 *
 * The controller is independent of the finite-element formulation. All model
 * operations are provided through callbacks, so the same class can operate on
 * constrained, reduced or otherwise transformed nonlinear systems.
 *
 * Stateful nonlinear subsystems are protected by increment trials that cover a
 * complete attempted arc-length increment and are committed only after the
 * increment is accepted.
 *
 * The optional active-set callback updates discontinuous nonlinear state after a
 * converged augmented Newton solve. If the active set changes, the controller
 * restarts Newton at the current load factor and path radius.
 */
class ArcLengthControl {
public:
    // Residual and tangent evaluation at the supplied reduced state and load
    // factor. The callback must overwrite both output objects.
    using Evaluate = std::function<void(
        const DynamicVector& q,
        Precision            lambda,
        DynamicVector&       residual,
        SparseMatrix&        tangent
    )>;

    // Linearized solve with a single right-hand side
    using LinearSolve = std::function<DynamicVector(
        const SparseMatrix&  tangent,
        const DynamicVector& rhs
    )>;

    // Linearized solve with one or more right-hand sides. Arc-length control
    // uses this for the augmented Newton system.
    using MatrixSolve = std::function<DynamicMatrix(
        const SparseMatrix&  tangent,
        const DynamicMatrix& rhs
    )>;

    // Problem-specific equilibrium-residual norm at the current load factor
    using ResidualNorm = std::function<Precision(
        const DynamicVector& residual,
        Precision            lambda
    )>;

    // Problem-specific correction norm for an accepted Newton correction
    using CorrectionNorm = std::function<Precision(
        const DynamicVector& q,
        const DynamicVector& dq
    )>;

    // Per-Newton-iteration reporting hook used by the owning load case
    using IterationCallback = std::function<void(
        Index     increment,
        Index     iteration,
        Precision lambda,
        Precision residual_norm,
        Precision correction_norm,
        Index     line_search_iterations,
        Time      assembly_ms,
        Time      solve_ms
    )>;

    // Accepted-increment reporting hook used for result writing and path storage
    using IncrementCallback = std::function<void(
        Index                increment,
        const DynamicVector& q,
        Precision            lambda
    )>;

    // Transaction callback for stateful nonlinear subsystems
    using TrialCallback = NewtonSolver::TrialCallback;

    // Optional post-convergence update for discontinuous nonlinear state. It
    // returns true when the active set is already consistent with the converged
    // state and false when Newton must be restarted on the same path increment.
    using ActiveSetCallback = std::function<bool(
        const DynamicVector& q,
        Precision            lambda
    )>;

public:
    ArcLengthControl() = default;

    bool solve(
        DynamicVector&           q,
        Precision&               lambda,
        const DynamicVector&     reference_load,
        const Evaluate&          evaluate,
        const LinearSolve&       linear_solve,
        const MatrixSolve&       matrix_solve,
        const ResidualNorm&      residual_norm,
        const CorrectionNorm&    correction_norm,
        const IterationCallback& on_iteration = {},
        const IncrementCallback& on_increment = {}
    );

    const char* failure_reason() const;

public:
    // Nonlinear increment and Newton limits
    Index maximum_increments = 100;
    Index maximum_iterations = 20;

    // Arc-length increment range, residual tolerance and load-factor weighting
    Precision tolerance         = Precision(1e-8);
    Precision initial_increment = Precision(0.1);
    Precision minimum_increment = Precision(1e-4);
    Precision maximum_increment = Precision(1.0);
    Precision psi               = Precision(1.0);

    // Adaptive increment controls
    Precision growth_factor     = Precision(1.5);
    Precision cutback_factor    = Precision(0.5);

    Index fast_iterations       = 6;
    Index slow_iterations       = 10;
    Index maximum_cutbacks      = 20;

    bool adaptive               = true;

    // Complete attempted-increment transaction
    TrialCallback begin_increment_trial;
    TrialCallback commit_increment_trial;
    TrialCallback rollback_increment_trial;

    // Optional discontinuous-state update after Newton convergence
    ActiveSetCallback update_active_set;

private:
    void reset_state_();
    void configure_newton_();
    void adapt_increment_();

    Precision arc_constraint_norm_(
        const DynamicVector& q,
        Precision            lambda
    ) const;

private:
    NewtonSolver newton_;

    Index     accepted_increments_ = 0;
    Precision increment_           = Precision(0);

    DynamicVector q_accepted_;
    Precision     lambda_accepted_ = Precision(0);

    DynamicVector previous_delta_q_;
    Precision     previous_delta_lambda_ = Precision(1);

    Precision load_scale2_  = Precision(0);
    Precision radius_scale_ = Precision(0);
    Precision arc_radius_   = Precision(0);

    bool adjusting_final_step_ = false;

    std::string failure_reason_ = "NONE";
};

} // namespace tools
} // namespace loadcase
} // namespace fem
