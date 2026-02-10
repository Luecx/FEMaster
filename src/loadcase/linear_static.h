/**
 * @file linear_static.h
 * @brief Declares the linear static load case.
 *
 * Solves the linear static equilibrium problem for the associated model using
 * configurable solver backends.
 *
 * @see src/loadcase/linear_static.cpp
 * @see src/solve/solver.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "loadcase.h"
#include "../solve/solver.h"

#include <string>
#include <vector>

namespace fem {
namespace loadcase {

/**
 * @struct LinearStatic
 * @brief Executes a linear static analysis on the model.
 */
struct LinearStatic : public LoadCase {
    std::vector<std::string> supps; ///< Support identifiers applied to the model.
    std::vector<std::string> loads; ///< Load identifiers applied to the model.
    solver::SolverDevice device = solver::CPU; ///< Solver device selection.
    solver::SolverMethod method = solver::DIRECT; ///< Solver method selection.
    std::string stiffness_file; ///< Optional path for stiffness matrix output.


    bool inertia_relief  = false; ///< Toggle for inertia relief (adds temporary inertial load to balance F/M).
	bool rebalance_loads = false; ///< Toggle for load rebalancing (adds loads so that sum F = sum M = 0).
    /**
     * @brief Constructs the linear static load case.
     *
     * @param id Load-case identifier.
     * @param writer Writer used for reporting.
     * @param model Model reference.
     */
    LinearStatic(ID id, reader::Writer* writer, model::Model* model);

    /**
     * @brief Executes the linear static solution procedure.
     */
    void run() override;
};

} // namespace loadcase
} // namespace fem
