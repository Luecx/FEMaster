/******************************************************************************
 * @file eigen_shift_invert_op_simple.cpp
 * @brief Implementation of y = (A - σ I)^{-1} x for the simple problem.
 *
 * @details Builds and factorizes the shifted matrix (A - σI) per call and then
 *          solves the linear system. MKL thread count is set from the global
 *          configuration when built with USE_MKL.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 ******************************************************************************/

#include "eigen_shift_invert_op_simple.h"
#include "../core/logging.h"

#include <cstring>   // std::memcpy

namespace {

// Ensure matrix is in compressed (CSC) format.
inline void ensure_compressed(const fem::SparseMatrix& M_const) {
    auto& M = const_cast<fem::SparseMatrix&>(M_const);
    if (!M.isCompressed()) M.makeCompressed();
}

} // unnamed namespace

namespace fem::solver {

ShiftInvertOpSimple::ShiftInvertOpSimple(const SparseMatrix& A, Scalar sigma)
    : _A_ref(A)
    , _sigma(sigma)
{
    ensure_compressed(_A_ref);
}

int ShiftInvertOpSimple::rows() const { return static_cast<int>(_A_ref.rows()); }
int ShiftInvertOpSimple::cols() const { return static_cast<int>(_A_ref.cols()); }

void ShiftInvertOpSimple::set_shift(const Scalar& s) {
    _sigma = s;
}

void ShiftInvertOpSimple::perform_op(const Scalar* x_in, Scalar* y_out) const {
#ifdef USE_MKL
    // Configure MKL threads from global config.
    mkl_set_num_threads(global_config.max_threads);
#endif

    const int n = static_cast<int>(_A_ref.rows());

    // Build M = A - σ I (in-place safe construction).
    SparseMatrix M = _A_ref;
    ensure_compressed(M);
    M.diagonal().array() -= _sigma;

    // Factorize and solve.
    _ldl.compute(M);
    if (_ldl.info() != Eigen::Success) {
        throw std::runtime_error("ShiftInvertOpSimple: LDLT factorization failed");
    }

    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x(x_in, n);
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y = _ldl.solve(x);

    std::memcpy(y_out, y.data(), sizeof(Scalar) * static_cast<size_t>(n));
}

} // namespace fem::solver
