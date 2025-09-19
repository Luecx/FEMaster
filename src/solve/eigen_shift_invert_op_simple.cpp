/******************************************************************************
 * @file eigen_shift_invert_op_simple.cpp
 * @brief Implementation of y = (A - σ I)^{-1} x with cached factorization.
 *
 * @details
 * Builds and caches the LDLT factorization of (A - σ I). If σ changes via
 * set_shift(), the next call to perform_op() will re-factorize.
 ******************************************************************************/

#include "eigen_shift_invert_op_simple.h"
#include "../core/logging.h"

#include <cstring>   // std::memcpy>
#include <stdexcept>

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
    // Defer re-factorization until next apply.
}

void ShiftInvertOpSimple::_factorize_if_needed() const
{
    if (std::isfinite(_sigma_fact) && _sigma_fact == _sigma) {
        return; // cached factorization is current
    }

    // Build M = A - σ I.
    SparseMatrix M = _A_ref;
    ensure_compressed(M);

    // Add (-σ) to the diagonal sparsely.
    const int n = static_cast<int>(M.rows());
    if (_sigma != Scalar(0)) {
        std::vector<Eigen::Triplet<Scalar>> trips;
        trips.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) trips.emplace_back(i, i, -_sigma);
        SparseMatrix D(n, n);
        D.setFromTriplets(trips.begin(), trips.end());
        M += D;
    }
    M.makeCompressed();

    // Factorize.
    _ldl.compute(M);
    if (_ldl.info() != Eigen::Success) {
        throw std::runtime_error("ShiftInvertOpSimple: LDLT factorization failed");
    }

    _sigma_fact = _sigma;
}

void ShiftInvertOpSimple::perform_op(const Scalar* x_in, Scalar* y_out) const
{
#ifdef USE_MKL
    // Configure MKL threads from global config each apply (cheap).
    mkl_set_num_threads(global_config.max_threads);
#endif

    _factorize_if_needed();

    const int n = static_cast<int>(_A_ref.rows());
    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x(x_in, n);
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y = _ldl.solve(x);
    if (_ldl.info() != Eigen::Success) {
        throw std::runtime_error("ShiftInvertOpSimple: solve() failed");
    }

    std::memcpy(y_out, y.data(), sizeof(Scalar) * static_cast<size_t>(n));
}

} // namespace fem::solver
