/******************************************************************************
* @file detect_rank.cpp
 * @brief Rank detection from diagonal of R with relative tolerance.
 ******************************************************************************/

#include "detect_rank.h"

#include "../../core/logging.h"
#include <cmath>

namespace fem { namespace constraint {

RankInfo detect_rank_from_R(const SparseMatrix& R, const RankSettings& s)
{
    RankInfo info{};

    const int rmax = std::min<int>(R.rows(), R.cols());
    Precision max_abs_diag = 0;

    for (int k = 0; k < rmax; ++k) {
        const Precision ad = std::abs(R.coeff(k, k));
        if (ad > max_abs_diag) max_abs_diag = ad;
    }
    info.max_abs_diag = max_abs_diag;

    const Precision tol = std::max<Precision>(Precision(1e-300), max_abs_diag * s.rel_tol);
    int r = 0;
    for (; r < rmax; ++r) {
        if (std::abs(R.coeff(r, r)) <= tol) break;
    }
    info.r = r;

    return info;
}

}} // namespace fem::constraint
