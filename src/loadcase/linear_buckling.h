/******************************************************************************
 * @file linear_buckling.h
 * @brief Linear buckling analysis using affine null-space constraints.
 *
 * @details
 * Solves the classical linear (eigenvalue) buckling problem:
 *
 *      K φ = λ (−K_g) φ,
 *
 * with constraints enforced via the affine null-space map u = u_p + T q.
 * Steps (high level):
 *   1) Pre-buckling static solve in reduced space (A q = b) to get displacements.
 *   2) From the pre-buckling field, compute stresses at integration points.
 *   3) Assemble geometric stiffness K_g (active DOFs).
 *   4) Reduce both K and K_g with T: A = Tᵀ K T, B = Tᵀ K_g T.
 *   5) Solve reduced generalized EVP: A φ = λ (−B) φ for the requested modes.
 *   6) Expand mode shapes back to node×6 layout for writing.
 *
 * Notes
 * -----
 * - This null-space route keeps the reduced A symmetric positive definite
 *   (assuming feasible constraints and SPD K), which is nice for solvers.
 * - Inhomogeneous constraints are handled via u_p in the preload solve;
 *   buckling itself is homogeneous (eigen problem).
 *
 * Example
 * -------
 * LinearBuckling buckling(1, writer, model, k=10);
 * buckling.supps = {"ENC_A", "TIE_1"};     // whatever your model expects
 * buckling.loads = {"PUSH_TOP"};           // preload to induce K_g
 * buckling.run();
 *
 * @author  Finn
 * @date    15.09.2025
 ******************************************************************************/

#pragma once

#include "loadcase.h"
#include "../solve/solver.h"

namespace fem { namespace loadcase {

struct LinearBuckling : public LoadCase {

    /**
     * @brief Construct a linear buckling load case.
     * @param id              Unique ID.
     * @param writer          Result writer.
     * @param model           FE model.
     * @param numEigenvalues  Number of modes to extract.
     */
    explicit LinearBuckling(ID id,
                            reader::Writer* writer,
                            model::Model* model,
                            int numEigenvalues);

    // User inputs
    std::vector<std::string> supps;           ///< Support/coupling identifiers → constraints.
    std::vector<std::string> loads;           ///< Load identifiers → preload (for K_g).
    int num_eigenvalues;                      ///< Number of buckling modes requested.
    Precision sigma = 0;                      ///< Target shift for eigenvalue search (0 = smallest).

    // === Debug / diagnostics ===
    /// If non-empty, write elastic stiffness (K, A) here.
    /// Will be suffixed "_K.mtx", "_A.mtx".
    std::string stiffness_file;

    /// If non-empty, write geometric stiffness (Kg, B) here.
    /// Will be suffixed "_Kg.mtx", "_B.mtx".
    std::string geom_file;

    // Solver selection
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;

    /// Execute the analysis.
    virtual void run() override;
};

}} // namespace fem::loadcase
