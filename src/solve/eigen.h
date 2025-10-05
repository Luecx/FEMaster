/**
 * @file eigen.h
 * @brief Sparse symmetric partial-eigen solvers (CPU, optional MKL/PARDISO).
 *
 * @details This header declares a compact API (backed by Spectra) to compute a
 *          small number (k) of extremal eigenpairs of large, **symmetric** and
 *          sparse problems:
 *
 *              • Simple:       A x = λ x
 *              • Generalized:  A x = λ B x
 *
 *          Supported modes:
 *              • Regular      – classic matrix–vector products (no solves)
 *              • ShiftInvert  – (A - σB)^{-1}B x   [simple: (A - σI)^{-1} x]
 *              • Buckling     – (A - σB)^{-1}A x   (returns physical λ)
 *              • Cayley       – (A - σB)^{-1}B x   (alternate transform)
 *
 *          Notes:
 *              • All paths assume symmetry; no unsymmetric/LU fallbacks.
 *              • With -DUSE_MKL, linear solves use MKL PardisoLDLT; otherwise
 *                Eigen's SimplicialLDLT is used.
 *              • GPU is not implemented here; requests fall back to CPU.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 */

#pragma once

#include "../core/core.h"     // Precision, DynamicVector/Matrix, SparseMatrix
#include "device.h"           // SolverDevice

#include <vector>
#include <utility>

namespace fem::solver {

/**
 * @struct EigenValueVectorPair
 * @brief A single eigenpair (λ, x).
 */
struct EigenValueVectorPair {
    Precision     value;   ///< Eigenvalue λ
    DynamicVector vector;  ///< Eigenvector x (size n)
};

/**
 * @enum EigenMode
 * @brief Selects the transformation/operator to be used by the solver.
 */
enum class EigenMode {
    Regular,      ///< No shift: standard operator (A or pair (A,B) with Cholesky)
    ShiftInvert,  ///< Simple: (A - σI)^{-1};  Generalized: (A - σB)^{-1}B
    Buckling,     ///< Generalized: (A - σB)^{-1}A   (physical λ)
    Cayley        ///< Generalized: (A - σB)^{-1}B   (alternate transform)
};

/**
 * @struct EigenOpts
 * @brief Execution parameters for the eigensolvers.
 */
struct EigenOpts {
    EigenMode mode = EigenMode::Regular; ///< Transform / operator
    double    sigma = 0.0;               ///< Shift (used by ShiftInvert/Buckling/Cayley)
    int       maxit = 3000;              ///< Max Arnoldi/Lanczos iterations
    double    tol   = 1e-8;              ///< Relative residual tolerance

    /// Sorting rule for returned eigenvalues.
    enum class Sort { LargestAlge, LargestMagn, SmallestAlge, SmallestMagn } sort = Sort::SmallestAlge;
};

/**
 * @brief Simple symmetric problem:  A x = λ x
 * @param device  CPU/GPU (GPU currently falls back to CPU)
 * @param A       Symmetric sparse matrix
 * @param k       Number of eigenpairs to compute (k > 0 and k < n)
 * @param opts    Algorithm options (mode must be Regular or ShiftInvert)
 * @return        Vector of k eigenpairs (sorted per `opts.sort`)
 */
std::vector<EigenValueVectorPair>
eigs(SolverDevice device,
     const SparseMatrix& A,
     int k,
     const EigenOpts& opts = {});

/**
 * @brief Generalized symmetric problem:  A x = λ B x
 * @param device  CPU/GPU (GPU currently falls back to CPU)
 * @param A       Symmetric sparse stiffness/left matrix
 * @param B       Symmetric sparse mass/right matrix
 * @param k       Number of eigenpairs to compute (k > 0 and k < n)
 * @param opts    Algorithm options (mode must be ShiftInvert, Buckling, or Cayley)
 * @return        Vector of k eigenpairs (sorted per `opts.sort`)
 */
std::vector<EigenValueVectorPair>
eigs(SolverDevice device,
     const SparseMatrix& A,
     const SparseMatrix& B,
     int k,
     const EigenOpts& opts);

} // namespace fem::solver
