/**
 * @file cb_compress.cpp
 * @brief Column compression of a sparse constraint matrix with optional row filtering.
 *
 * Algorithm (single pass over columns):
 *   1) Decide which columns survive: a column j survives if it has any entry
 *      in a row i that is kept (keep_row[i] == 1). If keep_row is empty, all rows are kept.
 *   2) Build `old2new` mapping and collect the surviving column indices in `used_cols`.
 *   3) Copy the surviving columns into `C_use`, skipping entries in dropped rows.
 *
 * Complexity:
 *   - O(nnz(C)) time, O(nnz(C_use)) memory for triplets.
 *   - `C_use` is returned compressed (makeCompressed()) and has the same number of rows as `C`.
 */

#include "./cb_compress.h"

namespace fem { namespace constraint {

void compress_zero_columns_with_row_filter(const SparseMatrix& C,
                                           const std::vector<char>& keep_row,
                                           std::vector<int>&   used_cols,
                                           Eigen::VectorXi&    old2new,
                                           SparseMatrix&       C_use)
{
    const int m = C.rows();
    const int n = C.cols();
    const bool has_filter = !keep_row.empty();

    // ---------------------------------------------------------------------
    // 1) Identify surviving columns (nonzero on kept rows)
    // ---------------------------------------------------------------------
    used_cols.clear();
    used_cols.reserve(n);

    for (int j = 0; j < C.outerSize(); ++j) {
        bool survives = false;

        if (!has_filter) {
            // Fast path: any nonzero means the column survives.
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                survives = true;
                break;
            }
        } else {
            // Filtered path: survive only if there exists an entry in a kept row.
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                if (keep_row[it.row()]) { survives = true; break; }
            }
        }

        if (survives) used_cols.push_back(j);
    }

    // Build remapping old -> new column index
    old2new = Eigen::VectorXi::Constant(n, -1);
    for (int k = 0; k < static_cast<int>(used_cols.size()); ++k) {
        old2new[ used_cols[k] ] = k;
    }

    // Early exit: if nothing survives, return an empty m×0 matrix.
    if (used_cols.empty()) {
        C_use.resize(m, 0);
        C_use.makeCompressed();
        return;
    }

    // ---------------------------------------------------------------------
    // 2) Copy surviving columns into C_use (skipping dropped rows)
    // ---------------------------------------------------------------------
    C_use.resize(m, static_cast<int>(used_cols.size()));

    // Heuristic reserve: we might over-reserve a bit; cheap and avoids reallocation.
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(static_cast<std::size_t>(C.nonZeros()));

    for (int jnew = 0; jnew < static_cast<int>(used_cols.size()); ++jnew) {
        const int j_old = used_cols[jnew];

        if (!has_filter) {
            // No row filter — copy the whole column
            for (SparseMatrix::InnerIterator it(C, j_old); it; ++it) {
                trips.emplace_back(it.row(), jnew, it.value());
            }
        } else {
            // Copy only entries located in kept rows
            for (SparseMatrix::InnerIterator it(C, j_old); it; ++it) {
                const int i = it.row();
                if (!keep_row[i]) continue;
                trips.emplace_back(i, jnew, it.value());
            }
        }
    }

    C_use.setFromTriplets(trips.begin(), trips.end());
    C_use.makeCompressed();
}

}} // namespace fem::constraint
