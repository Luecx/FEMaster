/******************************************************************************
 * @file preprocess.cpp
 * @brief Implements preprocessing of constraint systems prior to QR.
 *
 * Preprocessing identifies directly fixed DOFs, substitutes them into the
 * system, and compacts the remaining matrix rows and columns.
 *
 * @see src/constraints/builder/preprocess.h
 * @see src/constraints/builder/particular_solution.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "preprocess.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <algorithm>

namespace fem {
namespace constraint {
namespace {

/******************************************************************************
 * @brief Describes a single non-zero row encountered during scanning.
 ******************************************************************************/
struct SimpleRow {
    int row = -1;
    int col = -1;
    Precision coefficient = Precision(0);
    Precision rhs = Precision(0);
};

/******************************************************************************
 * @brief Finds rows with exactly one non-zero entry.
 ******************************************************************************/
std::vector<SimpleRow> find_single_nnz_rows(const SparseMatrix& C,
                                            const DynamicVector& d,
                                            bool homogeneous) {
    const int m = C.rows();
    const bool has_rhs = (d.size() == m) && !homogeneous;

    std::vector<int> nnz_row(m, 0);
    std::vector<int> last_col(m, -1);
    std::vector<Precision> last_val(m, Precision(0));

    for (int j = 0; j < C.outerSize(); ++j) {
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            const int i = it.row();
            nnz_row[i] += 1;
            last_col[i] = j;
            last_val[i] = it.value();
        }
    }

    std::vector<SimpleRow> singles;
    singles.reserve(m / 8 + 1);
    for (int i = 0; i < m; ++i) {
        if (nnz_row[i] == 1 && last_col[i] >= 0) {
            if (last_val[i] == Precision(0)) {
                continue;
            }
            SimpleRow sr;
            sr.row = i;
            sr.col = last_col[i];
            sr.coefficient = last_val[i];
            sr.rhs = has_rhs ? d[i] : Precision(0);
            singles.push_back(sr);
        }
    }
    return singles;
}

/******************************************************************************
 * @brief Compacts zero columns while respecting a row-keep mask.
 ******************************************************************************/
void compress_zero_columns_with_row_filter(const SparseMatrix& C,
                                           const std::vector<char>& keep_row,
                                           std::vector<int>& used_cols,
                                           Eigen::VectorXi& old2new,
                                           SparseMatrix& C_use) {
    const int n = C.cols();

    used_cols.clear();
    used_cols.reserve(n);

    for (int j = 0; j < C.outerSize(); ++j) {
        bool nonzero = false;
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            if (!keep_row.empty() && !keep_row[it.row()]) {
                continue;
            }
            nonzero = true;
            break;
        }
        if (nonzero) {
            used_cols.push_back(j);
        }
    }

    old2new = Eigen::VectorXi::Constant(n, -1);
    for (int k = 0; k < static_cast<int>(used_cols.size()); ++k) {
        old2new[used_cols[k]] = k;
    }

    C_use.resize(C.rows(), static_cast<int>(used_cols.size()));
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(C.nonZeros());

    for (int jnew = 0; jnew < static_cast<int>(used_cols.size()); ++jnew) {
        const int j = used_cols[jnew];
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            if (!keep_row.empty() && !keep_row[it.row()]) {
                continue;
            }
            trips.emplace_back(it.row(), jnew, it.value());
        }
    }
    C_use.setFromTriplets(trips.begin(), trips.end());
    C_use.makeCompressed();
}

} // namespace

/******************************************************************************
 * @copydoc preprocess_constraints
 ******************************************************************************/
PreprocessOutput preprocess_constraints(const PreprocessInput& input) {
    PreprocessOutput output;

    const int m = input.m;
    const int n = input.n;

    SparseMatrix C = input.C;
    DynamicVector d = input.d;

    auto singles = find_single_nnz_rows(C, d, input.homogeneous);

    output.is_fixed_col.assign(n, 0);
    output.fixed_val.assign(n, Precision(0));
    for (const auto& row : singles) {
        const Precision fixed_value = (row.coefficient != Precision(0)) ? (row.rhs / row.coefficient) : Precision(0);
        output.is_fixed_col[row.col] = 1;
        output.fixed_val[row.col] = fixed_value;
    }

    output.keep_row.assign(m, 1);
    for (const auto& row : singles) {
        output.keep_row[row.row] = 0;
    }

    if (!singles.empty()) {
        if (d.size() == 0) {
            d = DynamicVector::Zero(m);
        }

        for (int j = 0; j < n; ++j) {
            if (!output.is_fixed_col[j]) {
                continue;
            }
            const Precision fixed_value = output.fixed_val[j];
            if (fixed_value == Precision(0)) {
                continue;
            }
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                d[it.row()] -= it.value() * fixed_value;
            }
        }

        std::vector<Eigen::Triplet<Precision>> trips;
        trips.reserve(std::max<std::size_t>(1, C.nonZeros() - singles.size()));
        for (int j = 0; j < C.outerSize(); ++j) {
            if (output.is_fixed_col[j]) {
                continue;
            }
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                trips.emplace_back(it.row(), j, it.value());
            }
        }
        C.setZero();
        C.resize(m, n);
        C.setFromTriplets(trips.begin(), trips.end());
        C.makeCompressed();
    }

    compress_zero_columns_with_row_filter(C, output.keep_row, output.used, output.old2new, output.C_use);
    output.n_use = static_cast<int>(output.used.size());
    output.n_fixed = 0;
    for (char flag : output.is_fixed_col) {
        output.n_fixed += (flag ? 1 : 0);
    }

    output.d_mod = (d.size() == m) ? d : DynamicVector::Zero(m);

    return output;
}

} // namespace constraint
} // namespace fem
