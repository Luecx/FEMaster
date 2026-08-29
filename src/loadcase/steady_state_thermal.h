/**
 * @file steady_state_thermal.h
 * @brief Declares the linear steady-state thermal load case.
 *
 * The load case solves scalar heat conduction in thermally capable solid
 * elements with constant material conductivity. Prescribed temperatures are
 * selected from the model's common support collectors; heat sources, prescribed
 * fluxes and convection are intentionally outside the initial formulation.
 *
 * The resulting system is `K_T T = 0` subject to affine temperature constraints.
 * Conductivity assembly and heat-flux recovery remain responsibilities of
 * `model::Model` and `model::ThermalElement`.
 *
 * @see SteadyStateThermal
 * @see model::ThermalElement
 * @see bc::Temperature
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#pragma once

#include "loadcase.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../solve/sparse/solve_sparse.h"

#include <string>
#include <vector>

namespace fem::loadcase {

/**
 * @brief Solves linear stationary heat conduction with prescribed temperatures.
 *
 * Every node connected to a `model::ThermalElement` owns one scalar temperature
 * unknown. Element conductivity matrices form the symmetric global operator,
 * while `bc::Temperature` objects selected from the named support collectors
 * impose absolute nodal values through the common constraint transformer.
 *
 * No thermal load vector is assembled in this first formulation. A non-trivial
 * temperature field is driven entirely by differing prescribed boundary values.
 * Every disconnected thermal component must therefore contain at least one
 * temperature support to remove its constant-temperature null mode.
 *
 * The converged scalar nodal temperature and global integration-point heat flux
 * are written as `TEMPERATURE` and `HEAT_FLUX` result fields.
 */
struct SteadyStateThermal : LoadCase {
    // Common support collectors from which temperature definitions are selected
    std::vector<std::string> supps{};

    // Linear solver and affine constraint backends
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;
    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace;

    // Analysis identity and execution
    std::string type_name() const override { return "STEADYSTATETHERMAL"; }
    void run() override;
};

} // namespace fem::loadcase
