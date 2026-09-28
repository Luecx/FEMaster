#include "../src/loadcase/tools/newton_solver.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace fem;
using namespace fem::loadcase::tools;

namespace {

NewtonSolver::Evaluate constant_residual(Precision value) {
    return [value](const DynamicVector&,
                   DynamicVector& residual,
                   SparseMatrix& tangent) {
        residual = DynamicVector::Constant(1, value);
        tangent.resize(1, 1);
    };
}

NewtonSolver::Norm max_abs_norm() {
    return [](const DynamicVector& vector) {
        return vector.size() > 0
            ? vector.lpNorm<Eigen::Infinity>()
            : Precision(0);
    };
}

} // namespace

TEST(NewtonSolverConvergence, DoesNotAcceptOrdinaryForceToleranceOnFirstIteration) {
    NewtonSolver solver;
    solver.maximum_iterations   = 1;
    solver.residual_tolerance   = Precision(5e-3);
    solver.correction_tolerance = Precision(1e-2);
    solver.line_search_enabled  = false;
    solver.early_failure_detection = false;

    DynamicVector x = DynamicVector::Zero(1);
    Index linear_solves = 0;

    const bool converged = solver.solve(
        x,
        constant_residual(Precision(4e-3)),
        [&](const SparseMatrix&, const DynamicVector&) {
            ++linear_solves;
            return DynamicVector::Zero(1);
        },
        max_abs_norm(),
        [](const DynamicVector&, const DynamicVector&) {
            return Precision(0.2);
        }
    );

    EXPECT_FALSE(converged);
    EXPECT_EQ(linear_solves, 1);
    EXPECT_EQ(solver.iterations(), 1);
    EXPECT_STREQ(solver.failure_reason(), "MAXIMUM_ITERATIONS");
}

TEST(NewtonSolverConvergence, AcceptsPracticallyExactFirstResidualWithoutCorrection) {
    NewtonSolver solver;
    solver.maximum_iterations   = 1;
    solver.residual_tolerance   = Precision(5e-3);
    solver.correction_tolerance = Precision(1e-2);
    solver.line_search_enabled  = false;
    solver.early_failure_detection = false;

    DynamicVector x = DynamicVector::Zero(1);
    Index linear_solves = 0;

    const bool converged = solver.solve(
        x,
        constant_residual(Precision(1e-9)),
        [&](const SparseMatrix&, const DynamicVector&) {
            ++linear_solves;
            return DynamicVector::Zero(1);
        },
        max_abs_norm(),
        [](const DynamicVector&, const DynamicVector&) {
            return Precision(1);
        }
    );

    EXPECT_TRUE(converged);
    EXPECT_EQ(linear_solves, 0);
    EXPECT_EQ(solver.iterations(), 1);
}

TEST(NewtonSolverConvergence, RequiresSmallCorrectionAfterForceConvergence) {
    NewtonSolver solver;
    solver.maximum_iterations   = 3;
    solver.residual_tolerance   = Precision(5e-3);
    solver.correction_tolerance = Precision(1e-2);
    solver.line_search_enabled  = false;
    solver.early_failure_detection = false;

    DynamicVector x = DynamicVector::Zero(1);
    Index linear_solves = 0;

    const bool converged = solver.solve(
        x,
        [](const DynamicVector& state,
           DynamicVector& residual,
           SparseMatrix& tangent) {
            residual = DynamicVector::Constant(
                1,
                state(0) < Precision(0.5) ? Precision(1) : Precision(0)
            );
            tangent.resize(1, 1);
        },
        [&](const SparseMatrix&, const DynamicVector& residual) {
            ++linear_solves;
            return residual;
        },
        max_abs_norm(),
        [](const DynamicVector&, const DynamicVector& correction) {
            return correction.size() > 0
                ? correction.lpNorm<Eigen::Infinity>()
                : Precision(0);
        }
    );

    EXPECT_TRUE(converged);
    EXPECT_EQ(linear_solves, 2);
    EXPECT_EQ(solver.iterations(), 3);
    EXPECT_NEAR(x(0), Precision(1), Precision(1e-12));
}

TEST(NewtonSolverConvergence, FirstIterationShortcutRespectsStricterResidualTolerance) {
    NewtonSolver solver;
    solver.maximum_iterations   = 1;
    solver.residual_tolerance   = Precision(1e-10);
    solver.correction_tolerance = Precision(1e-2);
    solver.line_search_enabled  = false;
    solver.early_failure_detection = false;

    DynamicVector x = DynamicVector::Zero(1);
    Index linear_solves = 0;

    const bool converged = solver.solve(
        x,
        constant_residual(Precision(5e-9)),
        [&](const SparseMatrix&, const DynamicVector&) {
            ++linear_solves;
            return DynamicVector::Zero(1);
        },
        max_abs_norm(),
        [](const DynamicVector&, const DynamicVector&) {
            return Precision(0);
        }
    );

    EXPECT_FALSE(converged);
    EXPECT_EQ(linear_solves, 1);
}


TEST(NewtonSolverFailureDetection, CutsPersistentSlowConvergenceFromRateForecast) {
    NewtonSolver solver;
    solver.maximum_iterations              = 20;
    solver.residual_tolerance              = Precision(1e-3);
    solver.correction_tolerance            = Precision(1e-6);
    solver.line_search_enabled             = false;
    solver.early_failure_detection         = true;
    solver.convergence_check_start         = 4;
    solver.slow_convergence_check_start    = 8;
    solver.maximum_slow_convergence_checks = 3;

    const std::vector<Precision> residuals{
        Precision(1.636), Precision(1.356), Precision(1.045),
        Precision(0.9748), Precision(0.9210), Precision(0.8175),
        Precision(0.7135), Precision(0.6999), Precision(0.6994),
        Precision(0.6986), Precision(0.6977)
    };

    DynamicVector x = DynamicVector::Zero(1);
    Index evaluation = 0;

    const bool converged = solver.solve(
        x,
        [&](const DynamicVector&, DynamicVector& residual, SparseMatrix& tangent) {
            const Index index = std::min(
                evaluation++,
                static_cast<Index>(residuals.size() - 1)
            );
            residual = DynamicVector::Constant(1, residuals[index]);
            tangent.resize(1, 1);
        },
        [](const SparseMatrix&, const DynamicVector&) {
            return DynamicVector::Zero(1);
        },
        max_abs_norm(),
        [](const DynamicVector&, const DynamicVector&) {
            return Precision(1);
        }
    );

    EXPECT_FALSE(converged);
    EXPECT_EQ(solver.iterations(), 10);
    EXPECT_STREQ(solver.failure_reason(), "SLOW_CONVERGENCE");
}

TEST(NewtonSolverFailureDetection, AllowsLateAccelerationWithinIterationBudget) {
    NewtonSolver solver;
    solver.maximum_iterations              = 20;
    solver.residual_tolerance              = Precision(1e-3);
    solver.correction_tolerance            = Precision(1e-2);
    solver.line_search_enabled             = false;
    solver.early_failure_detection         = true;
    solver.convergence_check_start         = 4;
    solver.slow_convergence_check_start    = 8;
    solver.maximum_slow_convergence_checks = 3;

    const std::vector<Precision> residuals{
        Precision(1.418), Precision(1.238), Precision(1.158),
        Precision(1.093), Precision(0.9637), Precision(0.8530),
        Precision(0.7599), Precision(0.6721), Precision(0.4742),
        Precision(0.2669), Precision(0.06772), Precision(0.004413),
        Precision(2e-4)
    };

    DynamicVector x = DynamicVector::Zero(1);
    Index evaluation = 0;

    const bool converged = solver.solve(
        x,
        [&](const DynamicVector&, DynamicVector& residual, SparseMatrix& tangent) {
            const Index index = std::min(
                evaluation++,
                static_cast<Index>(residuals.size() - 1)
            );
            residual = DynamicVector::Constant(1, residuals[index]);
            tangent.resize(1, 1);
        },
        [](const SparseMatrix&, const DynamicVector&) {
            return DynamicVector::Zero(1);
        },
        max_abs_norm(),
        [](const DynamicVector&, const DynamicVector&) {
            return Precision(0);
        }
    );

    EXPECT_TRUE(converged);
    EXPECT_EQ(solver.iterations(), 13);
    EXPECT_STREQ(solver.failure_reason(), "NONE");
}

TEST(NewtonSolverFailureDetection, CutsRepeatedStrongLineSearchDampingWithoutProgress) {
    NewtonSolver solver;
    solver.maximum_iterations             = 20;
    solver.residual_tolerance             = Precision(1e-6);
    solver.correction_tolerance           = Precision(1e-6);
    solver.early_failure_detection        = true;
    solver.line_search_enabled            = true;
    solver.strong_damping_step_length     = Precision(1.0 / 64.0);
    solver.strong_damping_residual_ratio  = Precision(0.9);
    solver.maximum_strong_damping_steps   = 3;

    DynamicVector x = DynamicVector::Zero(1);
    Precision current_x = Precision(0);
    Precision current_residual = Precision(1);
    Index main_evaluations = 0;

    const bool converged = solver.solve(
        x,
        [&](const DynamicVector& state,
            DynamicVector& residual,
            SparseMatrix& tangent) {
            current_x = state(0);
            current_residual =
                std::pow(Precision(0.99), Precision(main_evaluations++));
            residual = DynamicVector::Constant(1, current_residual);
            tangent.resize(1, 1);
        },
        [](const SparseMatrix&, const DynamicVector&) {
            return DynamicVector::Ones(1);
        },
        max_abs_norm(),
        [](const DynamicVector&, const DynamicVector& correction) {
            return correction.lpNorm<Eigen::Infinity>();
        },
        {},
        [&](const DynamicVector& trial_state, DynamicVector& residual) {
            const Precision alpha = std::abs(trial_state(0) - current_x);
            const Precision value =
                alpha <= solver.strong_damping_step_length * Precision(1.0001)
                    ? Precision(0.99) * current_residual
                    : Precision(2) * current_residual;
            residual = DynamicVector::Constant(1, value);
        }
    );

    EXPECT_FALSE(converged);
    EXPECT_EQ(solver.iterations(), 4);
    EXPECT_STREQ(solver.failure_reason(), "LINE_SEARCH_STAGNATION");
}
