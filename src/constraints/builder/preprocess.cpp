/******************************************************************************
 * @file preprocess.cpp
 * @brief Implementation of sparse preprocessing: single-NNZ scan, substitution,
 *        and zero-column compression with row filter.
 ******************************************************************************/

#include "preprocess.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <algorithm>
#include <unordered_set>

namespace fem { namespace constraint {

namespace {

/** Single-NNZ row descriptor (internal). */
struct SimpleRow {
    int       row  = -1;        // row index in original C
    int       col  = -1;        // only nonzero column
    Precision a    = Precision(0);
    Precision d    = Precision(0);
};

/** Scan C to find rows with exactly one nonzero; read d if inhomogeneous. */
static std::vector<SimpleRow>
find_single_nnz_rows(const SparseMatrix& C, const DynamicVector& d, bool homogeneous)
{
    const int m = C.rows();
    const bool has_d = (d.size() == m) && !homogeneous;

    std::vector<int>      nnz_row(m, 0);
    std::vector<int>      last_col(m, -1);
    std::vector<Precision> last_val(m, Precision(0));

    for (int j = 0; j < C.outerSize(); ++j) {
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            const int i = it.row();
            nnz_row[i] += 1;
            last_col[i]  = j;
            last_val[i]  = it.value();
        }
    }

    std::vector<SimpleRow> singles;
    singles.reserve(m / 8 + 1);
    for (int i = 0; i < m; ++i) {
        if (nnz_row[i] == 1 && last_col[i] >= 0) {
            if (last_val[i] == Precision(0)) continue;
            SimpleRow sr;
            sr.row = i;
            sr.col = last_col[i];
            sr.a   = last_val[i];
            sr.d   = has_d ? d[i] : Precision(0);
            singles.push_back(sr);
        }
    }
    return singles;
}

/** Compress zero columns (considering a keep_row mask). */
static void compress_zero_columns_with_row_filter(const SparseMatrix& C,
                                                  const std::vector<char>& keep_row,
                                                  std::vector<int>&   used_cols,
                                                  Eigen::VectorXi&    old2new,
                                                  SparseMatrix&       C_use)
{
    const int m = C.rows();
    const int n = C.cols();

    used_cols.clear();
    used_cols.reserve(n);

    for (int j = 0; j < C.outerSize(); ++j) {
        bool nonzero = false;
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            if (!keep_row.empty() && !keep_row[it.row()]) continue;
            nonzero = true; break;
        }
        if (nonzero) used_cols.push_back(j);
    }

    old2new = Eigen::VectorXi::Constant(n, -1);
    for (int k = 0; k < (int)used_cols.size(); ++k) old2new[used_cols[k]] = k;

    C_use.resize(m, (int)used_cols.size());
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(C.nonZeros());

    for (int jnew = 0; jnew < (int)used_cols.size(); ++jnew) {
        const int j = used_cols[jnew];
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            if (!keep_row.empty() && !keep_row[it.row()]) continue;
            trips.emplace_back(it.row(), jnew, it.value());
        }
    }
    C_use.setFromTriplets(trips.begin(), trips.end());
    C_use.makeCompressed();
}

} // namespace


PreprocessOutput preprocess_constraints(const PreprocessInput& in)
{
    PreprocessOutput out{};

    const int m = in.m;
    const int n = in.n;

    SparseMatrix C = in.C;         // local mutable copy
    DynamicVector d = in.d;        // may be empty if homogeneous

    // --- 1) Single-NNZ scan --------------------------------------------------
    auto singles = find_single_nnz_rows(C, d, in.homogeneous);

    out.is_fixed_col.assign(n, 0);
    out.fixed_val.assign(n, Precision(0));
    for (const auto& sr : singles) {
        const Precision ufix = (sr.a != Precision(0)) ? (sr.d / sr.a) : Precision(0);
        out.is_fixed_col[sr.col] = 1;
        out.fixed_val[sr.col]    = ufix;
    }

    // Keep all rows except the single-NNZ equations (they are resolved already)
    out.keep_row.assign(m, 1);
    for (const auto& sr : singles) out.keep_row[sr.row] = 0;

    // --- 2) Substitute and zero fixed columns --------------------------------
    if (!singles.empty()) {
        if (d.size() == 0) d = DynamicVector::Zero(m);

        // d := d - C(:,j) * ufix
        for (int j = 0; j < n; ++j) {
            if (!out.is_fixed_col[j]) continue;
            const Precision ufix = out.fixed_val[j];
            if (ufix == Precision(0)) continue;
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                d[it.row()] -= it.value() * ufix;
            }
        }

        // Rebuild C skipping fixed columns (sparse-safe)
        std::vector<Eigen::Triplet<Precision>> trips;
        trips.reserve(std::max<size_t>(1, C.nonZeros() - singles.size()));
        for (int j = 0; j < C.outerSize(); ++j) {
            if (out.is_fixed_col[j]) continue;
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                trips.emplace_back(it.row(), j, it.value());
            }
        }
        C.setZero();
        C.resize(m, n);
        C.setFromTriplets(trips.begin(), trips.end());
        C.makeCompressed();
    }

    // --- 3) Compress columns with keep_row -----------------------------------
    compress_zero_columns_with_row_filter(C, out.keep_row, out.used, out.old2new, out.C_use);
    out.n_use   = (int)out.used.size();
    out.n_fixed = 0; for (char f : out.is_fixed_col) out.n_fixed += (f ? 1 : 0);

    // d_mod must be aligned to original row indexing (same size as m)
    out.d_mod = (d.size() == m) ? d : DynamicVector::Zero(m);

    return out;
}

}} // namespace fem::constraint
