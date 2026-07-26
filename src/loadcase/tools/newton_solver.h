/**
 * @file newton_solver.h
 * @brief Declares a generic Newton solver.
 */

#pragma once

#include "../../core/timer.h"
#include "../../core/types_eig.h"

#include <functional>

namespace fem {
namespace loadcase {
namespace tools {

class NewtonSolver {
public:
    using Evaluate = std::function<void(
        const DynamicVector& x,
        DynamicVector&       residual,
        SparseMatrix&        tangent
    )>;

    using LinearSolve = std::function<DynamicVector(
        const SparseMatrix&  tangent,
        const DynamicVector& residual
    )>;

    using Norm = std::function<Precision(
        const DynamicVector& vector
    )>;

    using CorrectionNorm = std::function<Precision(
        const DynamicVector& x,
        const DynamicVector& dx
    )>;

    using IterationCallback = std::function<void(
        Index     iteration,
        Precision residual_norm,
        Precision correction_norm,
        Precision convergence_order,
        Time      assembly_ms,
        Time      solve_ms,
        bool      converged
    )>;

    using TrialCallback = std::function<void()>;

public:
    NewtonSolver() = default;

    bool solve(
        DynamicVector&           x,
        const Evaluate&          evaluate,
        const LinearSolve&       linear_solve,
        const Norm&              residual_norm,
        const CorrectionNorm&    correction_norm,
        const IterationCallback& on_iteration = {}
    );

    Index     iterations()              const;
    Precision initial_residual_norm()   const;
    Precision last_residual_norm()      const;
    Precision last_correction_norm()    const;
    Precision convergence_order()       const;
    Precision last_step_length()        const;

    bool failed_by_divergence()         const;
    bool failed_by_residual_increase()  const;
    bool failed_by_stagnation()         const;
    bool failed_by_poor_reduction()     const;
    bool failed_by_line_search()        const;
    bool failed_by_maximum_iterations() const;

    const char* failure_reason() const;

public:
    Index maximum_iterations = 30;

    Precision residual_tolerance   = Precision(1e-8);
    Precision stagnation_tolerance = Precision(0);

    bool check_finite            = true;
    bool early_failure_detection = true;

    Index convergence_check_start = 4;

    Precision divergence_factor          = Precision(1e3);
    Precision minimum_residual_reduction = Precision(0.8);
    Precision stagnation_ratio           = Precision(1e-3);

    Index maximum_residual_increases = 2;
    Index maximum_stagnation_steps   = 2;

    bool      line_search_enabled             = true;
    Index     maximum_line_search_iterations  = 16;
    Precision line_search_reduction           = Precision(0.5);
    Precision line_search_sufficient_decrease = Precision(1e-4);
    Precision minimum_step_length             = Precision(1e-8);

    TrialCallback begin_line_search_trial;
    TrialCallback commit_line_search_trial;
    TrialCallback rollback_line_search_trial;

private:
    void reset_state_();

    void update_residual_history_();
    void update_convergence_order_();
    void update_failure_counters_();

    bool should_stop_early_();

private:
    Index iterations_ = 0;

    Precision initial_residual_norm_           = Precision(0);
    Precision last_residual_norm_              = Precision(0);
    Precision last_correction_norm_            = Precision(0);
    Precision previous_residual_norm_          = Precision(0);
    Precision previous_previous_residual_norm_ = Precision(0);
    Precision convergence_order_               = Precision(0);
    Precision last_step_length_                = Precision(1);

    Index residual_increase_count_ = 0;
    Index stagnation_count_        = 0;

    bool failed_by_divergence_         = false;
    bool failed_by_residual_increase_  = false;
    bool failed_by_stagnation_         = false;
    bool failed_by_poor_reduction_     = false;
    bool failed_by_line_search_        = false;
    bool failed_by_maximum_iterations_ = false;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
