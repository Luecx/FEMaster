/******************************************************************************
 * @file eigen_shift_invert_op_gen.cpp
 * @brief Implementation of y = (A - σ B)^{-1} x with cached factorization.
 *
 * @details
 * Builds and caches the LDLT factorization of (A - σ B). If σ changes via
 * set_shift(), the next call to perform_op() will re-factorize.
 ******************************************************************************/

#include "eigen_shift_invert_op_gen.h"
#include "../core/logging.h"

#include <cstring>    // std::memcpy
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
    // Defer re-factorization until next apply.
}

void ShiftInvertOpGeneral::_factorize_if_needed() const
{
    if (std::isfinite(_sigma_fact) && _sigma_fact == _sigma) {
        return; // cached factorization is current
    }

    // Build M = A - σ B (symmetric by assumption).
    SparseMatrix M = _A_ref;
    ensure_compressed(M);

    SparseMatrix SB = _B_ref;
    ensure_compressed(SB);

    if (_sigma != Scalar(0)) {
        M -= (_sigma * SB);
    }
    M.makeCompressed();

    // Factorize.
    _ldl.compute(M);
    if (_ldl.info() != Eigen::Success) {
        throw std::runtime_error("ShiftInvertOpGeneral: LDLT factorization failed");
    }
    logging::info(true, "Factorized (A - σ B) for σ = ", _sigma);

    _sigma_fact = _sigma;
}

void ShiftInvertOpGeneral::perform_op(const Scalar* x_in, Scalar* y_out) const
{
#ifdef USE_MKL
    // Configure MKL threads from global config each apply (cheap).
    mkl_set_num_threads(global_config.max_threads);
#endif

    // Basic dimension check to be defensive.
    if (_A_ref.rows() != _B_ref.rows() || _A_ref.cols() != _B_ref.cols()) {
        throw std::invalid_argument("ShiftInvertOpGeneral: A and B must have the same shape");
    }


    _factorize_if_needed();

    const int n = static_cast<int>(_A_ref.rows());
    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x(x_in, n);
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y = _ldl.solve(x);
    if (_ldl.info() != Eigen::Success) {
        throw std::runtime_error("ShiftInvertOpGeneral: solve() failed");
    }

    std::memcpy(y_out, y.data(), sizeof(Scalar) * static_cast<size_t>(n));
}

} // namespace fem::solver
