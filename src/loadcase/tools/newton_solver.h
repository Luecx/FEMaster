/**
 * @file newton_solver.h
 * @brief Defines a generic callback-based Newton solver.
 *
 * The Newton solver provides the nonlinear iteration kernel used by higher-level
 * solution strategies such as load control, displacement control and arc-length
 * methods.
 *
 * The solver itself is independent of the finite-element formulation. Residual
 * assembly, tangent assembly, linear-system solution and problem-specific norm
 * definitions are supplied through callbacks. This allows the same Newton
 * implementation to operate on full systems, reduced systems and transformed
 * systems without depending on a specific model or constraint representation.
 *
 * An optional backtracking line search can reduce the Newton step when the full
 * correction does not provide sufficient residual reduction. Trial-state
 * callbacks allow stateful nonlinear subsystems to begin, commit and roll back
 * temporary evaluations performed during the line search.
 *
 * The solver additionally records the diagnostics required for nonlinear path
 * control and can terminate early when it detects divergence, repeated residual
 * growth, residual stagnation or insufficient reduction from the initial
 * residual.
 *
 * Increment management, load-factor adaptation, active-set restarts and result
 * writing remain responsibilities of the calling solution strategy.
 *
 * @see NewtonSolver
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#pragma once

#include "../../core/timer.h"
#include "../../core/types_eig.h"

#include <functional>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * @brief Generic Newton solver for nonlinear algebraic systems.
 *
 * The solver seeks a state vector `x` satisfying a nonlinear equilibrium
 * equation
 *
 *     R(x) = 0.
 *
 * At every Newton iteration, the caller supplies the residual vector and
 * tangent matrix through the evaluation callback. The linear-solve callback
 * then computes the correction from the assembled linearized system.
 *
 * The sign convention of the Newton equation is deliberately not imposed by
 * this class. The solver updates the state according to
 *
 *     x <- x + dx,
 *
 * so the linear-solve callback must return the correction with the appropriate
 * sign for the residual convention used by the caller. For example, when
 *
 *     R = F_external - F_internal,
 *
 * the callback commonly solves
 *
 *     K_t dx = R.
 *
 * The residual and correction norms are also supplied by the caller. This is
 * important because nonlinear finite-element analyses often use scaled,
 * normalized or energy-based convergence measures rather than the raw
 * Euclidean vector norm.
 *
 * When line search is enabled, the solver evaluates trial states
 *
 *     x_trial = x + alpha dx,
 *
 * with progressively reduced step lengths `alpha`. The first trial satisfying
 * the sufficient-decrease condition is accepted. For non-smooth problems such
 * as contact, the best strictly residual-reducing trial may be accepted even
 * when the formal sufficient-decrease condition is not fulfilled.
 *
 * Trial callbacks provide transactional state handling for history-dependent
 * evaluations:
 *
 * - `begin_line_search_trial` prepares a temporary trial state,
 * - `commit_line_search_trial` accepts the evaluated trial state,
 * - `rollback_line_search_trial` restores the state before the trial.
 *
 * These callbacks are optional for purely state-independent problems but are
 * required when residual evaluation modifies material history or other
 * persistent trial data.
 *
 * A solver object may be reused for multiple nonlinear solves. Diagnostic and
 * failure state is reset at the beginning of every call to `solve`.
 */
class NewtonSolver {
public:
    // Residual and tangent evaluation at the supplied nonlinear state. The
    // callback must overwrite both output objects with the system associated
    // with the current state vector.
    using Evaluate = std::function<void(
        const DynamicVector& x,
        DynamicVector&       residual,
        SparseMatrix&        tangent
    )>;

    // Residual-only evaluation at a supplied nonlinear trial state. Line search
    // uses this optional callback when available because trial acceptance only
    // depends on the residual norm, not on a freshly assembled tangent.
    using EvaluateResidual = std::function<void(
        const DynamicVector& x,
        DynamicVector&       residual
    )>;

    // Linearized-system solution. The callback returns the correction applied
    // directly to the current state through x += dx.
    using LinearSolve = std::function<DynamicVector(
        const SparseMatrix&  tangent,
        const DynamicVector& residual
    )>;

    // Problem-specific residual convergence measure
    using Norm = std::function<Precision(
        const DynamicVector& vector
    )>;

    // Problem-specific correction convergence measure. The current state is
    // supplied in addition to the accepted correction so relative displacement
    // or solution norms can be evaluated.
    using CorrectionNorm = std::function<Precision(
        const DynamicVector& x,
        const DynamicVector& dx
    )>;

    // Per-iteration reporting callback. A converged iteration has no associated
    // linear solve and therefore reports a zero correction norm and zero solve
    // time.
    using IterationCallback = std::function<void(
        Index     iteration,
        Precision residual_norm,
        Precision correction_norm,
        Index     line_search_iterations,
        Time      assembly_ms,
        Time      solve_ms
    )>;

    // Transaction callback used to manage temporary state during line-search
    // residual evaluations
    using TrialCallback = std::function<void()>;

public:
    // Newton iteration limits and convergence tolerances
    Index     maximum_iterations   = 30;
    Precision residual_tolerance   = Precision(1e-8);
    Precision stagnation_tolerance = Precision(0);

    // Numerical validity and early-failure detection
    bool  check_finite            = true;
    bool  early_failure_detection = true;
    Index convergence_check_start = 4;

    // Early-failure thresholds. The minimum residual reduction is expressed as
    // the largest admissible ratio between the current and initial residual.
    Precision divergence_factor          = Precision(1e3);
    Precision minimum_residual_reduction = Precision(0.8);
    Precision stagnation_ratio           = Precision(1e-3);

    // Number of consecutive residual increases or stagnating residual changes
    // accepted before early termination
    Index maximum_residual_increases = 2;
    Index maximum_stagnation_steps   = 2;

    // Backtracking line-search settings
    bool      line_search_enabled             = true;
    Index     maximum_line_search_iterations  = 16;
    Precision line_search_reduction           = Precision(0.5);
    Precision line_search_sufficient_decrease = Precision(1e-4);
    Precision minimum_step_length             = Precision(1e-8);

    // Optional transactional state handling for line-search evaluations
    TrialCallback begin_line_search_trial;
    TrialCallback commit_line_search_trial;
    TrialCallback rollback_line_search_trial;

    // Construction
    NewtonSolver() = default;

    // Solve the nonlinear system starting from x. The state vector contains the
    // converged solution on success and the last accepted Newton state on
    // failure.
    bool solve(
        DynamicVector&           x,
        const Evaluate&          evaluate,
        const LinearSolve&       linear_solve,
        const Norm&              residual_norm,
        const CorrectionNorm&    correction_norm,
        const IterationCallback& on_iteration = {},
        const EvaluateResidual&  evaluate_residual = {}
    );

    // Diagnostics required by nonlinear path controllers
    Index     iterations()       const;
    Precision last_step_length() const;

    // Textual classification of the most recent failure
    const char* failure_reason() const;

private:
    // Iteration state of the most recent solve
    Index iterations_ = 0;

    Precision initial_residual_norm_  = Precision(0);
    Precision last_residual_norm_     = Precision(0);
    Precision last_correction_norm_   = Precision(0);
    Precision previous_residual_norm_ = Precision(0);
    Precision last_step_length_       = Precision(1);

    // Consecutive failure-detection counters
    Index residual_increase_count_ = 0;
    Index stagnation_count_        = 0;

    // Failure classification flags
    bool failed_by_divergence_         = false;
    bool failed_by_residual_increase_  = false;
    bool failed_by_stagnation_         = false;
    bool failed_by_poor_reduction_     = false;
    bool failed_by_line_search_        = false;
    bool failed_by_maximum_iterations_ = false;

    // Reset all diagnostics and failure state before a new solve
    void reset_state_();

    // Update residual history and failure-detection counters after each residual
    // evaluation
    void update_residual_history_();
    void update_failure_counters_();

    // Evaluate the configured early-failure criteria
    bool should_stop_early_();
};

} // namespace tools
} // namespace loadcase
} // namespace fem
