/**
 * @file solve_eigval_shiftinvert_op_simple.h
 * @brief Defines the Spectra shift-invert operator for standard symmetric eigenproblems.
 *
 * The operator applies
 *
 *     y = (A - sigma I)^-1 x
 *
 * during Spectra iterations. The shifted matrix is factorized lazily and the
 * factorization is cached while the shift remains unchanged, so repeated
 * operator applications only perform sparse triangular solves.
 *
 * The sparse factorization backend follows the FEMaster CPU configuration:
 * Intel oneMKL uses PardisoLDLT, Apple Accelerate uses AccelerateLDLT and the
 * portable fallback uses Eigen SimplicialLDLT. Spectra owns the iteration
 * strategy while this class is responsible only for the shifted linear solve.
 *
 * @author Finn Eggers
 * @date 19.09.2025
 */

#pragma once

#include "../../core/types_eig.h"

#include <Eigen/SparseCore>

#if defined(USE_MKL)
#include "../../core/config.h"
#include <Eigen/PardisoSupport>
#include <mkl.h>
#elif defined(USE_ACCELERATE)
#include <Eigen/AccelerateSupport>
#else
#include <Eigen/SparseCholesky>
#endif

#include <limits>

namespace fem::solver {

/**
 * @brief Spectra operator for standard symmetric shift-invert eigenvalue solves.
 *
 * The class stores a reference to the symmetric system matrix and the active
 * spectral shift. A factorization of `A - sigma I` is created on demand and
 * remains valid until `set_shift()` changes the requested shift.
 *
 * `perform_op()` therefore has the cost of a sparse solve during the regular
 * Spectra iteration and only triggers a new factorization after a shift change.
 * The concrete LDLT implementation is selected at compile time from oneMKL,
 * Apple Accelerate or Eigen.
 */
class ShiftInvertOpSimple {
public:
    using Scalar = Precision;

private:
    // Symmetric system matrix kept by reference for the complete Spectra solve
    const SparseMatrix& _A_ref;

    // Requested shift and the shift represented by the cached factorization
    Scalar         _sigma      = Scalar(0);
    mutable Scalar _sigma_fact = std::numeric_limits<Scalar>::quiet_NaN();

    // Cached factorization of A - sigma I
#if defined(USE_MKL)
    mutable Eigen::PardisoLDLT<SparseMatrix> _ldl;
#elif defined(USE_ACCELERATE)
    mutable Eigen::AccelerateLDLT<SparseMatrix> _ldl;
#else
    mutable Eigen::SimplicialLDLT<SparseMatrix> _ldl;
#endif

public:
    // Construction
    ShiftInvertOpSimple(const SparseMatrix& A, Scalar sigma);

    // Spectra operator dimensions
    int rows() const;
    int cols() const;

    // Shift control. Changing the shift invalidates the cached factorization
    // logically; the new factorization is built lazily on the next application.
    void set_shift(const Scalar& sigma);

    // Apply the inverse shifted operator to one Spectra vector
    void perform_op(const Scalar* x_in, Scalar* y_out) const;

private:
    // Refresh the cached factorization when it does not represent the active shift
    void _factorize_if_needed() const;
};

} // namespace fem::solver
