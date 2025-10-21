#pragma once
/**
 * @file
 * @brief Rayleigh damping model: C = α M + β K (minimal, header-only).
 *
 * Usage:
 *   #include "fem/damping/rayleigh.h"
 *   fem::damping::Rayleigh dmp(0.02, 1e-4);
 *   SparseMatrix C = dmp.build(M, K);          // C = αM + βK
 *   // If you need to accumulate:
 *   // C_total = C_total + dmp.build(M, K);
 *
 * Assumptions:
 *   - M, K are BC-reduced and constant (linear, time-invariant system).
 *   - Uses your project types from core/types_eig.h (SparseMatrix, Precision).
 */

#include "../core/types_eig.h"  // SparseMatrix, Precision

namespace fem {
namespace damping {

struct Rayleigh {
    Precision alpha = Precision(0);  ///< mass-proportional coefficient
    Precision beta  = Precision(0);  ///< stiffness-proportional coefficient

    Rayleigh() = default;
    Rayleigh(Precision a, Precision b) : alpha(a), beta(b) {}

    /// Build a fresh damping matrix: C = α M + β K
    [[nodiscard]] SparseMatrix
    build(const SparseMatrix& M, const SparseMatrix& K) const
    {
        // Fast paths to avoid unnecessary ops
        if (alpha == Precision(0) && beta == Precision(0)) {
            SparseMatrix C(M.rows(), M.cols());
            C.setZero();
            return C;
        }
        if (beta == Precision(0)) {
            return alpha * M;
        }
        if (alpha == Precision(0)) {
            return beta * K;
        }
        // General case
        SparseMatrix C = alpha * M;
        C = C + beta * K;
        return C;
    }
};

} // namespace damping
} // namespace fem
