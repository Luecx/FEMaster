/******************************************************************************
 * @file build_x.cpp
 * @brief Build sparse X = -R11^{-1} R12 via per-column upper-tri backsolves.
 ******************************************************************************/

#include "build_x.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <unordered_map>

namespace fem { namespace constraint {

XCols build_X_cols_from_R(const SparseMatrix& R, int r)
{
    XCols Xc{};

    const int ncols = R.cols();
    if (r <= 0 || ncols <= r) {
        Xc.cols.clear();
        return Xc;
    }

    const int nm_use = ncols - r;
    Xc.cols.assign(nm_use, {});

    // Build compact row-wise view of R11 (upper-triangular r×r)
    std::vector<Precision> diag(r, Precision(0));
    std::vector<std::vector<std::pair<int, Precision>>> row_upper(r);
    for (int j = 0; j < r; ++j) {
        for (SparseMatrix::InnerIterator it(R, j); it; ++it) {
            const int i = it.row();
            if (i >= r) continue;
            const Precision v = it.value();
            if (i == j) {
                diag[i] = v;
            } else if (i < j) {
                row_upper[i].push_back({j, v});
            }
        }
    }

    // Back-substitute each column j of R12 (global column = r + j)
    #pragma omp parallel for if (nm_use > 2)
    for (int j = 0; j < nm_use; ++j) {
        const int Rcol = r + j;

        // Build RHS: b = -R(0:r-1, Rcol)
        std::unordered_map<int, Precision> bmap;
        bmap.reserve(16);
        for (SparseMatrix::InnerIterator it(R, Rcol); it; ++it) {
            const int i = it.row();
            if (i >= r) continue;
            bmap[i] = -it.value();
        }
        if (bmap.empty()) continue;

        // Solve R11 x = b  (upper-triangular, back-substitution)
        std::vector<Precision> x(r, Precision(0));
        for (int i = r - 1; i >= 0; --i) {
            Precision sum = Precision(0);
            auto itb = bmap.find(i);
            if (itb != bmap.end()) sum += itb->second;

            const auto& row = row_upper[i];
            for (const auto& e : row) {
                const int k = e.first;
                sum -= e.second * x[k];
            }

            const Precision di = diag[i];
            x[i] = (di != Precision(0)) ? (sum / di) : Precision(0);
        }

        // Collect nonzeros for column j
        std::vector<std::pair<int, Precision>> col_entries;
        col_entries.reserve(16);
        for (int i = 0; i < r; ++i) {
            const Precision xi = x[i];
            if (xi != Precision(0)) col_entries.emplace_back(i, xi);
        }

        Xc.cols[j] = std::move(col_entries);
    }

    return Xc;
}

}} // namespace fem::constraint
