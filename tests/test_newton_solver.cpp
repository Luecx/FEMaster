#include "../src/loadcase/tools/newton_solver.h"

#include <gtest/gtest.h>

#include <cmath>

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
