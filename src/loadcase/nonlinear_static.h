/**
 * @file nonlinear_static.h
 * @brief Declares the nonlinear static load case.
 *
 * This file defines the nonlinear static load case used for incremental
 * equilibrium iterations. The implementation is based on an updated-Lagrangian
 * formulation and solves the nonlinear residual equation by repeated tangent
 * stiffness updates.
 */

#pragma once

#include "loadcase.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../solve/sparse/solve_sparse.h"

#include <string>
#include <vector>

namespace fem {
namespace loadcase {

/**
 * @brief Nonlinear path control method.
 */
enum class NonlinearControl {
    LoadControl,
    ArcLength
};

/**
 * @brief Nonlinear static load case using incremental equilibrium iterations.
 *
 * The nonlinear static load case solves a quasi-static problem by applying
 * loads over nonlinear increments. Within each increment, the equilibrium
 * residual is iteratively reduced until the configured convergence tolerance is
 * reached or the maximum number of iterations is exceeded.
 *
 * The current implementation is intended for updated-Lagrangian nonlinear
 * analyses, where the structural state is updated incrementally after each
 * converged load step.
 */
struct NonlinearStatic : public LoadCase {
    using ConstraintMethod = constraint::ConstraintTransformer::Method;

    /**
     * @brief Names of support collectors active in this load case.
     */
    std::vector<std::string> supps;

    /**
     * @brief Names of load collectors active in this load case.
     */
    std::vector<std::string> loads;

    /**
     * @brief Solver device used for the linearized systems.
     */
    solver::SolverDevice device = solver::CPU;

    /**
     * @brief Solver method used for the linearized systems.
     */
    solver::SolverMethod method = solver::DIRECT;

    /**
     * @brief Method used to enforce linear constraints.
     *
     * The default is a null-space transformation, i.e. the constrained problem
     * is projected into an unconstrained reduced coordinate space.
     */
    ConstraintMethod constraint_method = ConstraintMethod::NullSpace;

    /**
     * @brief Nonlinear path control method.
     *
     * LoadControl prescribes the load factor increment directly. ArcLength uses
     * the same increment value as an equivalent load-factor step and converts it
     * internally into an arc-length radius.
     */
    NonlinearControl control = NonlinearControl::LoadControl;

    /**
     * @brief Optional file path for exporting stiffness data.
     *
     * An empty string disables stiffness matrix output.
     */
    std::string stiffness_file;

    /**
     * @brief Maximum number of accepted nonlinear increments.
     *
     * For LoadControl this is mostly a safety limit because the analysis stops
     * once the load factor reaches one. For ArcLength it limits the number of
     * accepted path-following steps.
     */
    int max_increments = 100;

    /**
     * @brief Initial nonlinear step size.
     *
     * For LoadControl this is the initial load-factor increment. For ArcLength
     * this is an equivalent initial load-factor increment used to establish the
     * fixed arc-length radius scale from the first tangent predictor.
     */
    Precision initial_increment = Precision(0.1);

    /**
     * @brief Minimum allowed nonlinear step size.
     *
     * For LoadControl this limits the load-factor increment. For ArcLength this
     * limits the normalized arc-length radius relative to the fixed initial
     * predictor scale.
     */
    Precision minimum_increment = Precision(1e-4);

    /**
     * @brief Maximum allowed nonlinear step size.
     *
     * For LoadControl this limits the load-factor increment. For ArcLength this
     * limits the normalized arc-length radius relative to the fixed initial
     * predictor scale.
     */
    Precision maximum_increment = Precision(1.0);

    /**
     * @brief Enables automatic growth and cutback of nonlinear increments.
     *
     * When disabled, every attempted increment uses @ref initial_increment.
     */
    bool adaptive_increments = true;

    /**
     * @brief Multiplicative growth factor for quickly converged increments.
     */
    Precision growth_factor = Precision(1.5);

    /**
     * @brief Multiplicative cutback factor for rejected or slow increments.
     */
    Precision cutback_factor = Precision(0.5);

    /**
     * @brief Iteration count up to which an accepted increment is enlarged.
     */
    int fast_iterations = 6;

    /**
     * @brief Iteration count from which an accepted increment is reduced.
     */
    int slow_iterations = 10;

    /**
     * @brief Maximum number of consecutive cutbacks for one increment.
     */
    int maximum_cutbacks = 20;

    /**
     * @brief Maximum number of nonlinear iterations per increment.
     */
    int max_iterations = 20;

    /**
     * @brief Convergence tolerance for the nonlinear residual.
     *
     * The exact norm and normalization are defined by the implementation in
     * the corresponding source file.
     */
    Precision tolerance = Precision(1e-8);

    /**
     * @brief Weighting factor for the load-factor part of the arc-length constraint.
     *
     * The internal arc-length constraint is
     * ||Delta q||^2 + psi^2 * load_scale^2 * Delta lambda^2 = radius^2,
     * where load_scale = ||dq_load|| is fixed from the first predictor and
     * A dq_load = T^T f_total.
     */
    Precision arc_length_psi = Precision(1.0);

    /**
     * @brief Enables artificial regularization of zero-stiffness rows.
     *
     * This can improve robustness for temporarily unconstrained or mechanism-like
     * configurations, but should be used carefully because it modifies the
     * tangent system.
     */
    bool regularize_zero_stiffness_rows = true;

    /**
     * @brief Scaling factor for zero-stiffness row regularization.
     *
     * Only used when @ref regularize_zero_stiffness_rows is enabled.
     */
    Precision zero_stiffness_regularization_alpha = Precision(1e-4);

    /**
     * @brief Constructs a nonlinear static load case.
     *
     * @param id Load case id.
     * @param writer Result writer used for output fields.
     * @param model FEM model on which the load case operates.
     */
    NonlinearStatic(ID id,
                    reader::Writer* writer,
                    model::Model* model);

    /**
     * @brief Executes the nonlinear static analysis.
     */
    void run() override;
};

} // namespace loadcase
} // namespace fem
