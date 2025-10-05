/******************************************************************************
 * @file build_x.cpp
 * @brief Computes sparse columns of `X = -R11^{-1} R12` using back substitution.
 *
 * This routine extracts the blocks of the sparse QR factor required to build
 * the master-to-slave coupling matrix.
 *
 * @see src/constraints/builder/build_x.h
 * @see src/constraints/builder/detect_rank.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "build_x.h"

#include <omp.h>
#include <unordered_map>
#include <vector>

namespace fem {
namespace constraint {

/******************************************************************************
 * @copydoc build_X_cols_from_R
 ******************************************************************************/
XCols build_X_cols_from_R(const SparseMatrix& R, int r) {
    XCols cols;

    const int ncols = R.cols();
    if (r <= 0 || ncols <= r) {
        cols.cols.clear();
        return cols;
    }

    const int nm_use = ncols - r;
    cols.cols.assign(nm_use, {});

    std::vector<Precision> diag(r, Precision(0));
    std::vector<std::vector<std::pair<int, Precision>>> row_upper(r);
    for (int j = 0; j < r; ++j) {
        for (SparseMatrix::InnerIterator it(R, j); it; ++it) {
            const int i = it.row();
            if (i >= r) {
                continue;
            }
            const Precision value = it.value();
            if (i == j) {
                diag[i] = value;
            } else if (i < j) {
                row_upper[i].push_back({j, value});
            }
        }
    }

#pragma omp parallel for if (nm_use > 2)
    for (int j = 0; j < nm_use; ++j) {
        const int Rcol = r + j;
        std::unordered_map<int, Precision> bmap;
        bmap.reserve(16);
        for (SparseMatrix::InnerIterator it(R, Rcol); it; ++it) {
            const int i = it.row();
            if (i >= r) {
                continue;
            }
            bmap[i] = -it.value();
        }
        if (bmap.empty()) {
            continue;
        }

        std::vector<Precision> x(r, Precision(0));
        for (int i = r - 1; i >= 0; --i) {
            Precision sum = Precision(0);
            auto itb = bmap.find(i);
            if (itb != bmap.end()) {
                sum += itb->second;
            }

            const auto& row = row_upper[i];
            for (const auto& entry : row) {
                const int k = entry.first;
                sum -= entry.second * x[k];
            }

            const Precision di = diag[i];
            x[i] = (di != Precision(0)) ? (sum / di) : Precision(0);
        }

        std::vector<std::pair<int, Precision>> column_entries;
        column_entries.reserve(16);
        for (int i = 0; i < r; ++i) {
            const Precision xi = x[i];
            if (xi != Precision(0)) {
                column_entries.emplace_back(i, xi);
            }
        }

        cols.cols[j] = std::move(column_entries);
    }

    return cols;
}

} // namespace constraint
} // namespace fem
