/**
 * @file steady_state_heat_transfer.h
 * @brief Declares linear steady-state heat conduction with prescribed temperatures.
 *
 * The thermal primary variable is one scalar temperature per active node.
 * Material conductivity and linear surface convection form the stationary
 * operator; heat-flux and convection histories form its right-hand side.
 * Temperature boundary conditions are enforced by the common constraint
 * transformer. No thermal capacity or time integration is involved.
 *
 * @see model::ThermalElement
 * @see model::Model::build_thermal_conductivity_matrix
 */

#pragma once

#include "loadcase.h"
#include "../constraints/transformer/constraint_transformer.h"
#include "../solve/sparse/solve_sparse.h"

namespace fem::loadcase {

/**
 * @brief Solves the uncoupled stationary temperature field.
 *
 * The assembled equilibrium is (K_cond + K_film) T = q, constrained by
 * C T = d. All operators use the reference geometry; the result is an absolute
 * nodal temperature, not an increment. The converged field becomes the target
 * temperature committed by the common load-case lifecycle.
 */
struct SteadyStateHeatTransfer : LoadCase {
    using Method = constraint::ConstraintTransformer::Method;

    // Linear solver backend and scalar-temperature constraint strategy
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;
    Method constraint_method    = Method::Elimination;

    // End-of-step time used to evaluate thermal boundary amplitudes
    Precision step_period = Precision(1);


    // Analysis identity and execution
    std::string type_name() const override { return "STEADYHEATTRANSFER"; }
    void run() override;
};

} // namespace fem::loadcase
