/******************************************************************************
 * @file LinearEigenfrequency.h
 * @brief LinearEigenfrequency.h defines a load case for performing linear
 * eigenfrequency analysis. The class computes the natural frequencies and mode
 * shapes of a structure based on the finite element model under specified
 * boundary conditions.
 *
 * @details The analysis assumes linear elasticity and free vibration conditions.
 * It uses the finite element method to solve the eigenvalue problem resulting
 * from the mass and stiffness matrices of the model. The number of eigenvalues
 * (natural frequencies) to compute can be customized.
 *
 * @author Created by <Your Name>
 * all rights reserved
 * @date Created on <Creation Date>
 *
 ******************************************************************************/

#pragma once

#include "loadcase.h"
#include "../solve/solver.h"
#include "../math/interpolate.h"

namespace fem {
namespace loadcase {

/******************************************************************************
 * LinearEigenfrequency class
 * This class performs linear eigenfrequency analysis on a finite element model
 * to compute natural frequencies and mode shapes. It extends the LoadCase class,
 * implementing a run function for solving the eigenvalue problem.
 ******************************************************************************/
struct LinearEigenfrequency : public LoadCase {

    //-------------------------------------------------------------------------
    // Constructor
    //-------------------------------------------------------------------------
    /**
     * @brief Constructs a LinearEigenfrequency load case with a given ID,
     * writer, model, and the number of eigenvalues to compute.
     *
     * @param id The unique ID of the load case.
     * @param writer Pointer to a Writer object for output handling.
     * @param model Pointer to the finite element model.
     * @param numEigenvalues The number of eigenvalues (natural frequencies) to compute.
     */
    explicit LinearEigenfrequency(ID id, reader::Writer* writer, model::Model* model, int numEigenvalues);

    //-------------------------------------------------------------------------
    // Data Members
    //-------------------------------------------------------------------------
    std::vector<std::string> supps;  /**< List of support conditions applied to the model. */
    int num_eigenvalues; /**< Number of eigenvalues to compute in the analysis. */

public:

    //-------------------------------------------------------------------------
    // Run Analysis
    //-------------------------------------------------------------------------
    /**
     * @brief Executes the linear eigenfrequency analysis. This function overrides
     * the base run() method and computes the natural frequencies and mode shapes
     * of the model under the specified boundary conditions.
     */
    virtual void run() override;
};

} // namespace loadcase
} // namespace fem
