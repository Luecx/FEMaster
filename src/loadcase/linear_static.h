/******************************************************************************
 * @file LinearStatic.h
 * @brief LinearStatic.h defines a load case for performing linear static
 * analysis. The class computes displacements, stresses, and reactions for a
 * structure based on the finite element model under applied loads.
 *
 * @details The analysis assumes linear elasticity and static loading conditions.
 * It uses the finite element method to solve the system of equations resulting 
 * from the applied loads and boundary conditions. Solver options such as direct 
 * or indirect methods and device selection (CPU/GPU) can be customized.
 *
 * @author Created by <Your Name>
 * all rights reserved
 * @date Created on <Creation Date>
 *
 ******************************************************************************/

#pragma once

#include "loadcase.h"
#include "../solve/solver.h"

namespace fem {
namespace loadcase {

/******************************************************************************
 * LinearStatic class
 * This class performs linear static analysis on a finite element model to 
 * compute displacements, stresses, and reaction forces. It extends the LoadCase
 * class, implementing a run function for solving the static equilibrium equations.
 ******************************************************************************/
struct LinearStatic : public LoadCase {

    //-------------------------------------------------------------------------
    // Constructor
    //-------------------------------------------------------------------------
    /**
     * @brief Constructs a LinearStatic load case with a given ID, writer, and model.
     *
     * @param id The unique ID of the load case.
     * @param writer Pointer to a Writer object for output handling.
     * @param model Pointer to the finite element model.
     */
    explicit LinearStatic(ID id, reader::Writer* writer, model::Model* model);

    //-------------------------------------------------------------------------
    // Data Members
    //-------------------------------------------------------------------------
    std::vector<std::string> supps;  /**< List of support conditions applied to the model. */
    std::vector<std::string> loads;  /**< List of loads applied to the model. */

    solver::SolverDevice device = solver::CPU;    /**< Solver device used for computation (e.g., CPU, GPU). */
    solver::SolverMethod method = solver::DIRECT; /**< Solver method (e.g., DIRECT, INDIRECT). */

    std::string stiffness_file;
public:

    //-------------------------------------------------------------------------
    // Run Analysis
    //-------------------------------------------------------------------------
    /**
     * @brief Executes the linear static analysis. This function overrides
     * the base run() method and computes the displacements, stresses, and 
     * reaction forces in the model under the specified loads and boundary conditions.
     */
    virtual void run() override;
};

} // namespace loadcase
} // namespace fem
