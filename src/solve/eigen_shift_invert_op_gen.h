/******************************************************************************
 * @file eigen_shift_invert_op_gen.h
 * @brief Spectra-compatible operator: y = (A - σ B)^{-1} x (symmetric A,B).
 *
 * @details This operator is used by Spectra’s generalized shift–invert solver.
 *          On each `perform_op` call, it assembles the shifted matrix
 *          (A - σB), factorizes it with:
 *              • MKL PardisoLDLT if compiled with USE_MKL, or
 *              • Eigen SimplicialLDLT otherwise,
 *          and solves for y given x.
 *
 *          The factorization is intentionally not cached between calls
 *          (σ assumed to change).
 *
 * @note    A and B must be symmetric and of identical size.
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
  // Eigen backend: sparse LDLT
  #include <Eigen/SparseCholesky>
#else
  // MKL backend: PardisoLDLT
  #include <Eigen/PardisoSupport>
  #include "../core/config.h"   // global_config.max_threads
  #include <mkl.h>
#endif

namespace fem::solver {

/**
 * @class ShiftInvertOpGeneral
 * @brief Applies y = (A - σ B)^{-1} x for symmetric A,B.
 */
class ShiftInvertOpGeneral {
public:
    using Scalar = Precision;

    /** Construct with matrices A,B and initial shift σ. */
    ShiftInvertOpGeneral(const SparseMatrix& A,
                         const SparseMatrix& B,
                         Scalar sigma);

    /** Rows of the operator (equals A.rows()). */
    int rows() const;

    /** Cols of the operator (equals A.cols()). */
    int cols() const;

    /** Update the shift σ. Triggers re-factorization on next call. */
    void set_shift(const Scalar& sigma);

    /** Compute y = (A - σ B)^{-1} x. */
    void perform_op(const Scalar* x_in, Scalar* y_out) const;

private:
    // Left matrix A (symmetric), held by reference.
    const SparseMatrix& _A_ref;

    // Right matrix B (symmetric), held by reference.
    const SparseMatrix& _B_ref;

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
