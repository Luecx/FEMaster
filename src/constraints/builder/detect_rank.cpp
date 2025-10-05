/**
 * @file detect_rank.cpp
 * @brief Determines the numerical rank from a sparse QR factor.
 *
 * The implementation inspects diagonal entries of the `R` factor to identify
 * significant pivots relative to the specified tolerance.
 *
 * @see src/constraints/builder/detect_rank.h
 * @see src/constraints/builder/factorize_qr.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "detect_rank.h"

#include <algorithm>
#include <cmath>

namespace fem {
namespace constraint {

/**
 * @copydoc detect_rank_from_R
 */
RankInfo detect_rank_from_R(const SparseMatrix& R, const RankSettings& settings) {
    RankInfo info;

    const int rmax = std::min<int>(R.rows(), R.cols());
    Precision max_abs_diag = 0;
    for (int k = 0; k < rmax; ++k) {
        const Precision diag = std::abs(R.coeff(k, k));
        max_abs_diag = std::max(max_abs_diag, diag);
    }
    info.max_abs_diag = max_abs_diag;

    const Precision tol = std::max<Precision>(Precision(1e-300), max_abs_diag * settings.rel_tol);
    int rank = 0;
    for (; rank < rmax; ++rank) {
        if (std::abs(R.coeff(rank, rank)) <= tol) {
            break;
        }
    }
    info.r = rank;

    return info;
}

} // namespace constraint
} // namespace fem
