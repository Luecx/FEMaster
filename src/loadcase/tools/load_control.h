/**
 * @file load_control.h
 * @brief Declares a generic nonlinear load-control solver.
 */

#pragma once

#include "newton_solver.h"

#include <functional>

namespace fem {
namespace loadcase {
namespace tools {

class LoadControl {
public:
    using Evaluate = std::function<void(
        const DynamicVector& q,
        Precision            lambda,
        DynamicVector&       residual,
        SparseMatrix&        tangent
    )>;

    using LinearSolve = std::function<DynamicVector(
        const SparseMatrix&  tangent,
        const DynamicVector& rhs
    )>;

    using ResidualNorm = std::function<Precision(
        const DynamicVector& residual,
        Precision            lambda
    )>;

    using CorrectionNorm = std::function<Precision(
        const DynamicVector& q,
        const DynamicVector& dq
    )>;

    using IterationCallback = std::function<void(
        Index     increment,
        Index     iteration,
        Precision lambda,
        Precision residual_norm,
        Precision correction_norm,
        Precision convergence_order,
        Time      assembly_ms,
        Time      solve_ms,
        bool      converged
    )>;

    using IncrementCallback = std::function<void(
        Index                increment,
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
        const IncrementCallback& on_increment = {}
    );

    Index       accepted_increments() const;
    Precision   increment()           const;
    const char* failure_reason()      const;

public:
    Index maximum_increments = 100;
    Index maximum_iterations = 20;

    Precision tolerance         = Precision(1e-8);
    Precision initial_increment = Precision(0.1);
    Precision minimum_increment = Precision(1e-4);
    Precision maximum_increment = Precision(1.0);

    Precision growth_factor     = Precision(1.5);
    Precision cutback_factor    = Precision(0.5);

    Index fast_iterations       = 6;
    Index slow_iterations       = 10;
    Index maximum_cutbacks      = 20;

    bool adaptive               = true;

private:
    void reset_state_();
    void configure_newton_();
    void adapt_increment_();

private:
    NewtonSolver newton_;

    Index     accepted_increments_ = 0;
    Precision increment_           = Precision(0);

    const char* failure_reason_    = "NONE";
};

} // namespace tools
} // namespace loadcase
} // namespace fem
