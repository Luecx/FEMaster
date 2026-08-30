/**
 * @file steady_state_thermal.h
 * @brief Declares the linear steady-state thermal load case.
 *
 * The load case assembles scalar conductivity, prescribed heat flux and linear
 * convection over thermally capable elements. Prescribed temperatures are
 * selected from the common Dirichlet collectors and enforced through the common
 * constraint transformer.
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
 * One scalar temperature unknown is activated at every node of a thermal
 * element. Conductivity and convection form the system operator, while
 * prescribed heat flux and the ambient convection term form the right-hand
 * side. Selected temperature conditions define the affine Dirichlet system.
 */
struct SteadyStateThermal : LoadCase {
    // Selected temperature and thermal Neumann collectors
    std::vector<std::string> supps{};
    std::vector<std::string> loads{};

    // Linear solver and constraint enforcement settings
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;
    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace;

    // Analysis identity and execution
    std::string type_name() const override { return "STEADYSTATETHERMAL"; }
    void run() override;
};

} // namespace fem::loadcase
