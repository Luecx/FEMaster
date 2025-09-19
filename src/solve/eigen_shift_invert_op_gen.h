/******************************************************************************
 * @file eigen_shift_invert_op_gen.h
 * @brief Spectra-compatible operator: y = (A - σ B)^{-1} x for symmetric A,B.
 *
 * @details
 * This operator is used by Spectra’s generalized shift–invert solvers:
 * ShiftInvert, Buckling, and Cayley. It maintains a cached factorization of
 * (A - σ B) and only re-factorizes when σ changes.
 *
 * Backends:
 *   - If USE_MKL is defined: uses MKL PardisoLDLT.
 *   - Otherwise: uses Eigen SimplicialLDLT.
 *
 * Threading:
 *   - When built with MKL, the MKL thread count is set from
 *     global_config.max_threads inside perform_op().
 *
 * Usage pattern with Spectra:
 *   - Spectra calls perform_op() many times with constant σ per solve, so
 *     caching the factorization is critical for performance.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date Created on 19.09.2025
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include <Eigen/SparseCore>

#ifndef USE_MKL
  #include <Eigen/SparseCholesky>     // SimplicialLDLT
#else
  #include <Eigen/PardisoSupport>     // PardisoLDLT
  #include "../core/config.h"         // global_config.max_threads
  #include <mkl.h>
#endif

#include <limits>   // std::numeric_limits

namespace fem::solver {

/******************************************************************************
 * @class ShiftInvertOpGeneral
 * @brief Applies y = (A - σ B)^{-1} x with cached factorization.
 ******************************************************************************/
class ShiftInvertOpGeneral {
public:
    using Scalar = Precision;

    /**
     * @brief Construct with matrices A,B and initial shift σ.
     * @param A     Reference to symmetric matrix A (kept by reference).
     * @param B     Reference to symmetric matrix B (kept by reference).
     * @param sigma Initial shift value.
     */
    ShiftInvertOpGeneral(const SparseMatrix& A,
                         const SparseMatrix& B,
                         Scalar sigma);

    /** @return number of rows (equals A.rows()). */
    int rows() const;

    /** @return number of columns (equals A.cols()). */
    int cols() const;

    /**
     * @brief Update the shift σ. The next apply will re-factorize.
     * @param sigma New shift value.
     */
    void set_shift(const Scalar& sigma);

    /**
     * @brief Compute y = (A - σ B)^{-1} x.
     * @param x_in  Pointer to input vector x (size = rows()).
     * @param y_out Pointer to output vector y (size = rows()).
     */
    void perform_op(const Scalar* x_in, Scalar* y_out) const;

private:
    // Matrix A (by reference).
    const SparseMatrix& _A_ref;

    // Matrix B (by reference).
    const SparseMatrix& _B_ref;

    // Target shift σ requested by the solver.
    Scalar _sigma{0};

    // Shift value that the cached factorization currently corresponds to.
    // Initialized to NaN to force the first factorization.
    mutable Scalar _sigma_fact{std::numeric_limits<Scalar>::quiet_NaN()};

#ifdef USE_MKL
    // MKL backend factorization object (LDLT of (A - σ B)).
    mutable Eigen::PardisoLDLT<SparseMatrix> _ldl;
#else
    // Eigen backend factorization object (SimplicialLDLT of (A - σ B)).
    mutable Eigen::SimplicialLDLT<SparseMatrix> _ldl;
#endif

    // Ensure a valid factorization exists for the current _sigma.
    void _factorize_if_needed() const;
};

} // namespace fem::solver
