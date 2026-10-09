/**
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
 */

#pragma once

#include "loadcase.h"
#include "../solve/sparse/solve_sparse.h"

namespace fem {
namespace loadcase {

/**
 * LinearEigenfrequency class
 * This class performs linear eigenfrequency analysis on a finite element model
 * to compute natural frequencies and mode shapes. It extends the LoadCase class,
 * implementing a run function for solving the eigenvalue problem.
 *
 * The reduced eigenproblem (T^T K T) phi = lambda (T^T M T) phi uses the
 * supports and model kinematic constraints currently active when run() begins.
 * Eigenvalue count/range, solver choices and mode output belong to the analysis;
 * condition history and condition activation belong to ModelData. External load
 * history remains stored in the model and is not used as an eigenproblem RHS.
 */
struct LinearEigenfrequency : public LoadCase {
    //-------------------------------------------------------------------------
    // Data Members
    //-------------------------------------------------------------------------
    int num_eigenvalues = 10; /**< Number of eigenvalues to compute in the analysis. */
    bool      use_eigenvalue_range = false;
    Precision min_eigenvalue       = 0;
    Precision max_eigenvalue       = 0;

    // Solver selection
    solver::SolverDevice device = solver::CPU;    ///< CPU / GPU.
    solver::SolverMethod method = solver::DIRECT; ///< DIRECT / INDIRECT - always DIRECT.

    // Construction and default result requests
    LinearEigenfrequency() {
        using io::writer::OutputField;
        output.set_defaults({
            OutputField::MODE_SHAPE
        });
    }

public:
    // Analysis identity and execution
    std::string type_name() const override { return "EIGENFREQ"; }
    void run() override;
};
} // namespace loadcase
} // namespace fem
