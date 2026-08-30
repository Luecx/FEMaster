/**
 * @file steady_state_thermal.h
 * @brief Declares the linear steady-state thermal load case.
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
 * @brief Solves stationary heat conduction with Dirichlet and Neumann data.
 */
struct SteadyStateThermal : LoadCase {
    // Dirichlet support collectors and thermal Neumann load collectors.
    std::vector<std::string> supps{};
    std::vector<std::string> loads{};

    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;
    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace;

    std::string type_name() const override { return "STEADYSTATETHERMAL"; }
    void run() override;
};

} // namespace fem::loadcase
