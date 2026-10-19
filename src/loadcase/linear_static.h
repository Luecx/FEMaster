
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
 * Solves K_T(u0,T0) du = (f_ext,end-f_ext,start) - (f_int(u0,T_target)-f_int(u0,T0))
 * about the prepared base with the configured constraint transformer and solver.
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
    using Method = constraint::ConstraintTransformer::Method;

    // Solver backend and constraint enforcement method
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;
    Method constraint_method    = Method::Elimination;

    // Optional file path for exporting the assembled stiffness matrices
    std::string stiffness_file;

    // Physical step duration used for evaluating prescribed load amplitudes
    Precision step_period = Precision(1);

    // Stores the actual external nodal load at the end of the last committed
    // static analysis. Linear perturbations use this field as the reference
    // load to determine the incremental loading (f_target - f_start).
    // Unlike the mutable load definitions, this field preserves the complete
    // assembled load, including previously activated collectors and the final
    // load factor of a nonlinear analysis. A null pointer represents zero
    // reference loading.
    model::Field::Ptr step_start_load = nullptr;

    // Optional free-body equilibrium methods. Inertia relief introduces
    // equivalent inertial forces to balance the applied loads and can include
    // concentrated point masses. Load rebalancing instead adjusts the applied
    // loads to eliminate the resultant force and moment.
    bool inertia_relief                       = false;
    bool inertia_relief_consider_point_masses = true;
    bool rebalance_loads                      = false;

    // Initialize the default structural output requests
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
