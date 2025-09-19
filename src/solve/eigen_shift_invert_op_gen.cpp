/******************************************************************************
 * @file eigen_shift_invert_op_gen.cpp
 * @brief Implementation of y = (A - σ B)^{-1} x for the generalized problem.
 *
 * @details Builds and factorizes the shifted matrix (A - σB) per call and then
 *          solves the linear system. MKL thread count is set from the global
 *          configuration when built with USE_MKL. A and B are required to be
 *          symmetric and of identical shape.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 ******************************************************************************/

#include "eigen_shift_invert_op_gen.h"
#include "../core/logging.h"

#include <cstring>   // std::memcpy
#include <stdexcept>

namespace {

// Ensure matrix is in compressed (CSC) format.
inline void ensure_compressed(const fem::SparseMatrix& M_const) {
    auto& M = const_cast<fem::SparseMatrix&>(M_const);
    if (!M.isCompressed()) M.makeCompressed();
}

} // unnamed namespace

namespace fem::solver {

ShiftInvertOpGeneral::ShiftInvertOpGeneral(const SparseMatrix& A,
                                           const SparseMatrix& B,
                                           Scalar sigma)
    : _A_ref(A)
    , _B_ref(B)
    , _sigma(sigma)
{
    ensure_compressed(_A_ref);
    ensure_compressed(_B_ref);
}

int ShiftInvertOpGeneral::rows() const { return static_cast<int>(_A_ref.rows()); }
int ShiftInvertOpGeneral::cols() const { return static_cast<int>(_A_ref.cols()); }

void ShiftInvertOpGeneral::set_shift(const Scalar& s) {
    _sigma = s;
}

void ShiftInvertOpGeneral::perform_op(const Scalar* x_in, Scalar* y_out) const {
#ifdef USE_MKL
    // Configure MKL threads from global config.
    mkl_set_num_threads(global_config.max_threads);
#endif

    const int nA = static_cast<int>(_A_ref.rows());
    const int nB = static_cast<int>(_B_ref.rows());
    if (nA != nB || _A_ref.cols() != _B_ref.cols()) {
        throw std::invalid_argument("ShiftInvertOpGeneral: A and B must have the same shape");
    }
    const int n = nA;

    // Build M = A - σ B.
    SparseMatrix M = _A_ref;
    ensure_compressed(M);
    SparseMatrix SB = _B_ref;
    ensure_compressed(SB);
    M -= (_sigma * SB);
    M.makeCompressed();

    // Factorize and solve.
    _ldl.compute(M);
    if (_ldl.info() != Eigen::Success) {
        throw std::runtime_error("ShiftInvertOpGeneral: LDLT factorization failed");
    }

    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x(x_in, n);
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y = _ldl.solve(x);

    std::memcpy(y_out, y.data(), sizeof(Scalar) * static_cast<size_t>(n));
}

} // namespace fem::solver
