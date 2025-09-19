/******************************************************************************
 * @file eigen_shift_invert_op_simple.h
 * @brief Spectra-compatible operator: y = (A - σ I)^{-1} x (symmetric A).
 *
 * @details This operator is used by Spectra’s shift–invert solver for the
 *          simple eigenproblem. On each `perform_op` call, it assembles the
 *          shifted matrix (A - σI), factorizes it with:
 *              • MKL PardisoLDLT (when compiled with USE_MKL), or
 *              • Eigen SimplicialLDLT (fallback),
 *          and solves for y given x.
 *
 *          The factorization is intentionally not cached between calls because
 *          σ is assumed to change across iterations.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
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

namespace fem::solver {

class ShiftInvertOpSimple {
public:
    using Scalar = Precision;

    /** Construct with matrix A and initial shift σ. */
    ShiftInvertOpSimple(const SparseMatrix& A, Scalar sigma);

    /** Rows of the operator (equals A.rows()). */
    int rows() const;

    /** Cols of the operator (equals A.cols()). */
    int cols() const;

    /** Update the shift σ. Triggers re-factorization on next call. */
    void set_shift(const Scalar& sigma);

    /** Compute y = (A - σ I)^{-1} x. */
    void perform_op(const Scalar* x_in, Scalar* y_out) const;

private:
    // Original symmetric matrix A (by reference).
    const SparseMatrix& _A_ref;

    // Current shift σ.
    Scalar _sigma{0};

#ifdef USE_MKL
    // MKL backend factorization object.
    mutable Eigen::PardisoLDLT<SparseMatrix> _ldl;
#else
    // Eigen fallback factorization object.
    mutable Eigen::SimplicialLDLT<SparseMatrix> _ldl;
#endif
};

} // namespace fem::solver
