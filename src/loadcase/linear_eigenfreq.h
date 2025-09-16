/******************************************************************************
* @file LinearEigenfrequency.h
 * @brief Linear eigenfrequency analysis using affine null-space constraints.
 *
 * Solves the generalized EVP
 *     K u = λ M u
 * with constraints enforced via the null-space map u = u_p + T q,
 * yielding the reduced EVP
 *     (Tᵀ K T) φ = λ (Tᵀ M T) φ .
 *
 * Outputs eigenvalues λ, natural frequencies f = √λ / (2π), mode shapes,
 * and simple modal participation factors in the 6 global DOF directions.
 *
 * @date    15.09.2025
 * @author  Finn
 ******************************************************************************/

#pragma once

#include "loadcase.h"
#include "../solve/solver.h"

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

    // Solver selection
    solver::SolverDevice device = solver::CPU;    ///< CPU / GPU.
    solver::SolverMethod method = solver::DIRECT; ///< DIRECT / INDIRECT - always DIRECT.
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
