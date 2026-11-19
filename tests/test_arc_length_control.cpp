/**
 * @file test_arc_length_control.cpp
 * @brief Regression tests for configuration-dependent Riks load derivatives.
 */

#include "../src/loadcase/tools/arc_length_control.h"

#include <Eigen/LU>
#include <gtest/gtest.h>
#include <cmath>

using namespace fem;

TEST(ArcLengthControl, CurrentLoadDerivativeIsUsedByPredictorAndBorderedNewton) {
    // One-DOF equilibrium with a configuration-dependent external load:
    // r(q,lambda) = lambda*(1 + alpha*q) - q.
    // The exact load-factor derivative is 1 + alpha*q, while the effective
    // stiffness is 1 - alpha*lambda. At lambda=1 the solution is q=1/(1-alpha).
    constexpr Precision alpha = Precision(0.2);

    loadcase::tools::ArcLengthControl controller;
    controller.psi = Precision(1);
    controller.tolerance = Precision(1e-9);
    controller.maximum_iterations = 30;

    auto& increments = controller.increment_manager;
    increments.initial_increment = Precision(0.05);
    increments.minimum_increment = Precision(1e-6);
    increments.maximum_increment = Precision(0.15);
    increments.maximum_increments = 100;
    increments.maximum_cutbacks = 30;

    Index derivative_calls = 0;
    Index predictor_solves = 0;
    Index bordered_solves = 0;
    bool saw_deformed_load_direction = false;
    Precision last_derivative = Precision(0);

    controller.load_derivative =
        [&](const DynamicVector& q, Precision) -> DynamicVector {
            ++derivative_calls;
            last_derivative = Precision(1) + alpha * q(0);
            saw_deformed_load_direction |= std::abs(q(0)) > Precision(1e-5);

            DynamicVector result(1);
            result(0) = last_derivative;
            return result;
        };

    auto evaluate = [&](const DynamicVector& q, Precision lambda,
                        DynamicVector& residual, SparseMatrix& tangent) {
        residual.resize(1);
        residual(0) = lambda * (Precision(1) + alpha * q(0)) - q(0);

        tangent.resize(1, 1);
        tangent.setZero();
        tangent.coeffRef(0, 0) = Precision(1) - alpha * lambda;
    };

    auto linear_solve = [&](const SparseMatrix& tangent,
                            const DynamicVector& rhs) -> DynamicVector {
        ++predictor_solves;
        // The controller must use the most recently evaluated derivative
        // for its tangent predictor, not the fixed reference load.
        EXPECT_NEAR(rhs(0), last_derivative, 1e-12);
        DynamicMatrix dense(tangent);
        return dense.fullPivLu().solve(rhs);
    };

    auto matrix_solve = [&](const SparseMatrix& tangent,
                            const DynamicMatrix& rhs) -> DynamicMatrix {
        ++bordered_solves;
        // The upper-right column of the bordered matrix is -df_ext/dlambda.
        EXPECT_NEAR(tangent.coeff(0, 1), -last_derivative, 1e-12);
        DynamicMatrix dense(tangent);
        return dense.fullPivLu().solve(rhs);
    };

    const auto residual_norm = [](const DynamicVector& residual, Precision) {
        return residual.lpNorm<Eigen::Infinity>();
    };
    const auto correction_norm = [](const DynamicVector&, const DynamicVector& dq) {
        return dq.norm();
    };

    DynamicVector q = DynamicVector::Zero(1);
    DynamicVector reference_load = DynamicVector::Ones(1);
    Precision lambda = Precision(0);

    ASSERT_TRUE(controller.solve(
        q, lambda, reference_load, evaluate, linear_solve, matrix_solve,
        residual_norm, correction_norm));
    EXPECT_NEAR(lambda, Precision(1), 1e-8);
    EXPECT_NEAR(q(0), Precision(1) / (Precision(1) - alpha), 1e-7);
    EXPECT_GT(predictor_solves, 1);
    EXPECT_GT(bordered_solves, 0);
    EXPECT_GT(derivative_calls, bordered_solves);
    EXPECT_TRUE(saw_deformed_load_direction);
}
