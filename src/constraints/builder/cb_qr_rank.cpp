#include "cb_qr_rank.h"
#include <algorithm>
#include <cmath>

namespace fem { namespace constraint {

void factor_and_rank_inplace(const SparseMatrix& C_use,
                             Precision rank_tol_rel,
                             QRRank& out)
{
    // Configure column pivot threshold (Eigen interprets this relatively, [0..1])
    const Precision thr = (rank_tol_rel <= 0 ? 0 : std::min<Precision>(0.1, rank_tol_rel));
    out.sqr.setPivotThreshold(thr);
    out.sqr.compute(C_use);
    out.R = out.sqr.matrixR();

    // Rank detection using relative diagonal test
    const int rmax = std::min<int>(out.R.rows(), out.R.cols());
    Precision max_abs_diag = 0;
    for (int k = 0; k < rmax; ++k) {
        max_abs_diag = std::max(max_abs_diag, std::abs(out.R.coeff(k, k)));
    }
    out.R11_max_diag = max_abs_diag;

    const Precision tol = std::max<Precision>(Precision(1e-300), max_abs_diag * rank_tol_rel);
    int r = 0;
    for (; r < rmax; ++r) {
        if (std::abs(out.R.coeff(r, r)) <= tol) break;
    }
    out.r = r;
}

}} // namespace fem::constraint
