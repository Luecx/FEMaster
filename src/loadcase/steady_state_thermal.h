/**
 * @file steady_state_thermal.h
 * @brief Declares the linear steady-state thermal load case.
 *
 * The load case assembles scalar element conductivity together with natural and
 * mixed thermal boundary conditions over thermally capable model topology.
 * Prescribed temperatures are selected from common Dirichlet collectors and
 * enforced through the common constraint transformer. Thermal load collectors
 * may contain pure Neumann heat-flow conditions and Robin conditions such as
 * convection.
 *
 * @see SteadyStateThermal
 * @see model::ThermalElement
 * @see bc::Temperature
 * @see bc::HeatFlux
 * @see bc::Convection
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "loadcase.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../solve/sparse/solve_sparse.h"

#include <string>
#include <vector>

namespace fem::loadcase {

/**
 * @brief Solves stationary linear heat conduction with thermal boundary data.
 *
 * One scalar temperature unknown is activated at every node connected to a
 * thermal element. Element conductivity and Robin boundary operators form the
 * left-hand side, while prescribed Neumann heat flow and the source part of
 * Robin conditions form the right-hand side. Selected `Temperature` conditions
 * define affine Dirichlet equations in the same scalar nodal unknown space.
 *
 * The loadcase reuses the generic constraint transformer and sparse linear-solver
 * infrastructure, but keeps its temperature DOF mapping separate from the
 * six-component structural system mapping.
 */
struct SteadyStateThermal : LoadCase {
    // Selected Dirichlet support collectors containing prescribed temperatures.
    // Incompatible structural support entries are filtered by the model's typed
    // thermal constraint collection path.
    std::vector<std::string> supps{};

    // Selected thermal load collectors. Compatible entries may be pure thermal
    // Neumann conditions or Robin conditions contributing to both RHS and LHS.
    std::vector<std::string> loads{};

    // Linear-system backend and constraint-enforcement method used for the scalar
    // steady-state thermal equilibrium problem.
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;
    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace;

    // Loadcase identity exposed to parser/diagnostic infrastructure and the main
    // execution entry point implemented in `steady_state_thermal.cpp`.
    std::string type_name() const override { return "STEADYSTATETHERMAL"; }
    void run() override;
};

} // namespace fem::loadcase
