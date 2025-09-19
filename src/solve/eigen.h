/******************************************************************************
 * @file eigen.h
 * @brief Sparse symmetric partial-eigen solvers (CPU, optional MKL/PARDISO).
 *
 * This header declares a tiny API around the Spectra library to compute
 * a handful (k) of extremal eigenpairs of large, **symmetric** sparse problems:
 *
 *   • Simple:       A x = λ x
 *   • Generalized:  A x = λ B x
 *
 * Modes:
 *   • Regular      – classic matrix–vector products (no linear solves)
 *   • ShiftInvert  – y = (A - σB)^{-1}(B x)  [simple: (A - σI)^{-1} x]
 *   • Buckling     – y = (A - σB)^{-1}(A x)  (returns physical λ)
 *   • Cayley       – y = (A - σB)^{-1}(B x)  (alternate transform)
 *
 * Notes
 *   • All paths assume symmetry. No unsymmetric fallbacks are provided.
 *   • If compiled with -DUSE_MKL, shift–invert uses MKL PARDISO (LDLT).
 *   • GPU is not implemented here; requests fall back to CPU with a note.
 *
 * @author
 *   Clean readability pass
 * @date 19.09.2025
 ******************************************************************************/

#pragma once

#include "../core/core.h"     // brings in Precision, DynamicVector/Matrix, SparseMatrix
#include "device.h"           // SolverDevice

#include <vector>
#include <utility>

namespace fem::solver {

// -----------------------------------------------------------------------------
// Return type
// -----------------------------------------------------------------------------
struct EigenValueVectorPair {
    Precision     value;   ///< Eigenvalue λ
    DynamicVector vector;  ///< Eigenvector x (size n)
};

// -----------------------------------------------------------------------------
// Controls (kept minimal; ncv is always auto-chosen internally)
// -----------------------------------------------------------------------------
enum class EigenMode {
    Regular,      ///< No shift
    ShiftInvert,  ///< (A - σB)^{-1}B (or (A-σI)^{-1} for simple)
    Buckling,     ///< (A - σB)^{-1}A   (requires A symmetric; typically A ≻ 0)
    Cayley        ///< Alternate transform; same inner solve as ShiftInvert
};

struct EigenOpts {
    EigenMode mode = EigenMode::Regular;
    double    sigma = 0.0;        ///< Used for ShiftInvert/Buckling/Cayley
    int       maxit = 3000;
    double    tol   = 1e-8;
    enum class Sort { LargestAlge, LargestMagn, SmallestAlge, SmallestMagn } sort = Sort::SmallestAlge;
};

// -----------------------------------------------------------------------------
// Simple problem:  A x = λ x  (A symmetric)
// -----------------------------------------------------------------------------
std::vector<EigenValueVectorPair>
eigs(SolverDevice device,
     const SparseMatrix& A,
     int k,
     const EigenOpts& opts = {});

// -----------------------------------------------------------------------------
// Generalized problem:  A x = λ B x  (A,B symmetric; B typically SPD)
// -----------------------------------------------------------------------------
std::vector<EigenValueVectorPair>
eigs(SolverDevice device,
     const SparseMatrix& A,
     const SparseMatrix& B,
     int k,
     const EigenOpts& opts);

} // namespace fem::solver
