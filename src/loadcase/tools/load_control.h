/**
 * @file load_control.h
 * @brief Declares direct load-factor control for nonlinear equilibrium paths.
 *
 * Load control advances a nonlinear static analysis by prescribed increments of
 * the scalar load factor. Within every accepted load increment, a generic Newton
 * solver reduces the residual at the fixed target load factor.
 *
 * The class owns only increment scheduling, cutback/growth decisions, active-set
 * restarts and callback wiring. Finite-element assembly, constraint projection,
 * sparse solves, convergence norms, result writing and stateful nonlinear
 * subsystem management are supplied by the calling load case.
 *
 * Trial callbacks are intentionally generic. They are used by contact today and
 * can also protect later history-dependent material states, damage variables or
 * other nonlinear state that must be committed only after an accepted step.
 *
 * @see LoadControl
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
 * @brief Increment controller for fixed-load-factor nonlinear solves.
 *
 * `LoadControl` attempts to move the scalar load factor from its current value
 * to one by repeatedly solving
 *
 *     R(q, lambda_target) = 0
 *
 * in the reduced nonlinear unknowns `q`. The next target value is obtained from
 * the current load factor plus the current increment size. Successful increments
 * may grow or shrink the following increment according to the configured
 * iteration-count thresholds; failed increments are rolled back and retried with
 * a cutback when adaptive stepping is enabled.
 *
 * The controller does not know the finite-element model, the residual
 * definition, the tangent matrix, the linear constraints or the sparse solver.
 * All those operations are provided through callbacks. This keeps the class
 * usable for transformed systems and for different nonlinear formulations.
 *
 * Stateful nonlinear subsystems are handled through two nested transaction
 * levels:
 *
 * - increment trials cover a complete attempted load increment and are committed
 *   only after the increment is accepted,
 * - line-search trials cover temporary Newton step-length evaluations and are
 *   committed only when the trial step is accepted by the Newton solver.
 *
 * The optional active-set callback updates discontinuous state after a converged
 * Newton solve. If the callback reports a changed active set, Newton is restarted
 * at the same target load factor using the newly committed state.
 */
class LoadControl {
public:
    // Residual and tangent evaluation at the supplied reduced state and load
    // factor. The callback must overwrite both output objects.
    using Evaluate = std::function<void(
        const DynamicVector& q,
        Precision            lambda,
        DynamicVector&       residual,
        SparseMatrix&        tangent
    )>;

    // Residual-only evaluation used by Newton line search when the tangent is
    // not needed for trial-step acceptance.
    using EvaluateResidual = std::function<void(
        const DynamicVector& q,
        Precision            lambda,
        DynamicVector&       residual
    )>;

    // Linearized Newton solve. The returned vector is applied directly as
    // q += dq by the generic Newton solver.
    using LinearSolve = std::function<DynamicVector(
        const SparseMatrix&  tangent,
        const DynamicVector& rhs
    )>;

    // Problem-specific residual norm at the current target load factor
    using ResidualNorm = std::function<Precision(
        const DynamicVector& residual,
        Precision            lambda
    )>;

    // Problem-specific correction norm for an accepted Newton correction
    using CorrectionNorm = std::function<Precision(
        const DynamicVector& q,
        const DynamicVector& dq
    )>;

    // Optional predictor that initializes the nonlinear state for a new load
    // increment before Newton iterations begin.
    using Predictor = std::function<void(
        DynamicVector& q,
        Precision      lambda,
        Precision      target_lambda
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
    // state and false when Newton must be restarted at the same load factor.
    using ActiveSetCallback = std::function<bool(
        const DynamicVector& q,
        Precision            lambda
    )>;

public:
    LoadControl() = default;

    bool solve(
        DynamicVector&           q,
        Precision&               lambda,
        const Evaluate&          evaluate,
        const LinearSolve&       linear_solve,
        const ResidualNorm&      residual_norm,
        const CorrectionNorm&    correction_norm,
        const IterationCallback& on_iteration = {},
        const IncrementCallback& on_increment = {},
        const Predictor&         predictor = {},
        const EvaluateResidual&  evaluate_residual = {}
    );

    const char* failure_reason() const;

public:
    // Nonlinear increment and Newton limits
    Index maximum_increments = 100;
    Index maximum_iterations = 20;

    // Load-factor increment range and residual tolerance
    Precision tolerance         = Precision(1e-8);
    Precision initial_increment = Precision(0.1);
    Precision minimum_increment = Precision(1e-4);
    Precision maximum_increment = Precision(1.0);

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

    // Temporary line-search transaction nested inside the current increment
    TrialCallback begin_line_search_trial;
    TrialCallback commit_line_search_trial;
    TrialCallback rollback_line_search_trial;

    // Optional discontinuous-state update after Newton convergence
    ActiveSetCallback update_active_set;

private:
    void reset_state_();
    void configure_newton_();
    void adapt_increment_();

private:
    NewtonSolver newton_;

    Index     accepted_increments_ = 0;
    Precision increment_           = Precision(0);

    std::string failure_reason_    = "NONE";
};

} // namespace tools
} // namespace loadcase
} // namespace fem
