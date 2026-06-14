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
 * @brief Nonlinear static load case using incremental equilibrium iterations.
 *
 * The nonlinear static load case solves a quasi-static problem by applying
 * loads over a prescribed number of increments. Within each increment, the
 * equilibrium residual is iteratively reduced until the configured convergence
 * tolerance is reached or the maximum number of iterations is exceeded.
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
     * @brief Optional file path for exporting or importing stiffness data.
     *
     * The precise interpretation depends on the nonlinear static implementation.
     * An empty string disables stiffness file output/input.
     */
    std::string stiffness_file;

    /**
     * @brief Number of load increments.
     *
     * The total external load is applied incrementally. Larger values generally
     * improve robustness for strongly nonlinear problems, but increase runtime.
     */
    int num_increments = 10;

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