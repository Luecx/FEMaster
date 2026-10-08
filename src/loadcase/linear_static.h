/**
 * @file linear_static.h
 * @brief Declares the linear static load case.
 *
 * Solves the linear static equilibrium problem for the associated model using
 * configurable solver backends.
 *
 * @see src/loadcase/linear_static.cpp
 * @see src/solve/sparse/solve_sparse.h
 * @author Finn Eggers
 * @date 06.03.2025
 *
 * Solves K u = f with the configured constraint transformer and sparse backend.
 * Loads and prescribed displacements are assembled from the current ModelData
 * history and collector-derived activations during run(). Only solver, output and free-body
 * balancing settings belong to this analysis; reusable condition definitions
 * remain model-owned. Temporary inertia-relief RBMs are removed after constraint
 * collection and do not modify persistent support history.
 */

#pragma once

#include "loadcase.h"
#include "../constraints/transformer/constraint_transformer.h"
#include "../solve/sparse/solve_sparse.h"

#include <string>

namespace fem {
namespace loadcase {

/**
 * @struct LinearStatic
 * @brief Executes a linear static analysis on the model.
 */
struct LinearStatic : public LoadCase {
    solver::SolverDevice device = solver::CPU; ///< Solver device selection.
    solver::SolverMethod method = solver::DIRECT; ///< Solver method selection.
    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace; ///< Constraint backend selection.
    std::string stiffness_file; ///< Optional path for stiffness matrix output.
    Precision step_period = Precision(1); ///< Step-end time for evaluating prescribed load amplitudes.
    bool inertia_relief  = false; ///< Toggle for inertia relief (adds temporary inertial load to balance F/M).
    bool inertia_relief_consider_point_masses = true; ///< Include POINTMASS features in inertia-relief mass/inertia and load assembly.
	bool rebalance_loads = false; ///< Toggle for load rebalancing (adds loads so that sum F = sum M = 0).

    // Construction and default result requests
    LinearStatic() {
        using io::writer::OutputField;
        output.set_defaults({
            OutputField::DISPLACEMENT,
            OutputField::STRAIN,
            OutputField::STRESS,
            OutputField::STRESS_TOP,
            OutputField::STRESS_BOT,
            OutputField::SHELL_RESULTANTS,
            OutputField::EXTERNAL_FORCES,
            OutputField::REACTION_FORCES,
            OutputField::LOCAL_SECTION_FORCES,
            OutputField::SHEAR_FLOW
        });
    }


    // Analysis identity and execution
    std::string type_name() const override { return "LINEARSTATIC"; }
    void run() override;
};
} // namespace loadcase
} // namespace fem
