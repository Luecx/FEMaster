/**
 * @file cb_build_x_sparse.cpp
 * @brief Sparse construction of Xd = -R11^{-1} R12 from SparseQR's R, without dense blocks.
 *
 * Overview
 * --------
 * We are given the upper-triangular factor R from a rank-revealing SparseQR on C_use * P = Q * R.
 * After partitioning columns by rank r:
 *   R = [ R11  R12  R13? ]
 *       [  0    0    0  ]   (rows below r are irrelevant here)
 *
 * We need Xd = -R11^{-1} R12, but must avoid forming R11/R12 as dense blocks.
 * This implementation:
 *   (1) Builds a compact row-major view of R11 (diagonal + strictly-upper row entries).
 *   (2) For each column j of R12 (j = 0..nm_use-1), forms the sparse RHS b = -R12(:,j),
 *       and back-substitutes through R11 to get the sparse column Xd(:,j).
 *
 * Notes
 * -----
 * - R is accessed by column via Eigen::InnerIterator; we ignore rows >= r.
 * - Each Xd column is solved independently and written into out.X_cols[j] (thread-safe).
 * - No dense r×nm_use temporaries are created.
 */

#include "cb_build_x_sparse.h"
#include <unordered_map>   // sparse RHS per column
#include <algorithm>       // std::max

namespace fem { namespace constraint {

XBuildOut build_X_sparse(const SparseMatrix& Rsp, int r, int nm_use)
{
    XBuildOut out;
    out.R11_diag.assign(r, Precision(0));
    out.R11_rows.assign(r, {});
    out.X_cols.assign(nm_use, {});

    // Trivial cases: nothing to do
    if (r == 0 || nm_use == 0) return out;

    // -------------------------------------------------------------------------
    // 1) Build a compact row-wise view of R11 (top-left r×r block of R).
    //
    // For each column j < r, iterate its nonzeros and keep only rows i < r.
    //  - Diagonal: store in R11_diag[i]
    //  - Strictly-upper entries (i < j): append to R11_rows[i] as (k=j, a_ik)
    //
    // This produces, per row i:
    //   R11(i, i) in R11_diag[i]
    //   { (k, R11(i,k)) for k=i+1..r-1 } in R11_rows[i]
    // -------------------------------------------------------------------------
    const int Rcols = Rsp.cols();
    for (int j = 0; j < Rcols && j < r; ++j) {
        for (SparseMatrix::InnerIterator it(Rsp, j); it; ++it) {
            const int i = it.row();
            if (i >= r) continue;           // below the R11 block
            const Precision v = it.value();
            if (i == j) {
                out.R11_diag[i] = v;
            } else if (i < j) {
                out.R11_rows[i].push_back({ j, v });
            } // i > j cannot happen (R is upper-triangular)
        }
    }

    // Optional: we could sort each R11_rows[i] by k to improve cache behavior.
    // Not strictly required since back-sub loops exactly over these entries.

    // -------------------------------------------------------------------------
    // 2) For each column j in R12, solve R11 * x = b where b = -R12(:,j).
    //
    // Implementation details:
    //  - Collect b sparsely by scanning column (r + j) of R and keeping rows < r.
    //  - Back-substitute from i = r-1 down to 0:
    //       x[i] = ( b[i] - sum_{k>i} R11(i,k)*x[k] ) / R11(i,i)
    //    If a diagonal entry is numerically zero, we defensively treat x[i]=0.
    //  - Finally, compress x to a sparse list of nonzeros (i, x[i]).
    //
    // Threading:
    //  - Each j writes only to out.X_cols[j], so the loop can be parallelized.
    // -------------------------------------------------------------------------
    #pragma omp parallel for if (nm_use > 2)
    for (int j = 0; j < nm_use; ++j) {
        const int Rcol = r + j;  // column in R that corresponds to R12[:, j]

        // Build sparse RHS b = -R12(:, j)
        std::unordered_map<int, Precision> bmap;
        bmap.reserve(16);
        for (SparseMatrix::InnerIterator it(Rsp, Rcol); it; ++it) {
            const int i = it.row();
            if (i >= r) continue;       // only rows in the R11 block matter
            bmap[i] = -it.value();
        }
        if (bmap.empty()) {
            // Column is zero ⇒ Xd(:,j) = 0
            continue;
        }

        // Dense temporary for the back-substitution
        std::vector<Precision> x(r, Precision(0));

        // Backward solve through the upper-triangular R11
        for (int i = r - 1; i >= 0; --i) {
            Precision acc = Precision(0);

            // Start with b[i] if present
            if (auto itb = bmap.find(i); itb != bmap.end()) {
                acc += itb->second;
            }

            // Subtract Σ_{k>i} R11(i,k) * x[k]
            const auto& row = out.R11_rows[i];
            for (const auto& e : row) {
                acc -= e.a * x[e.k];
            }

            const Precision di = out.R11_diag[i];
            x[i] = (di != Precision(0)) ? (acc / di) : Precision(0);
        }

        // Compress to sparse entries (skip exact zeros)
        std::vector<std::pair<int, Precision>> col_entries;
        col_entries.reserve(16);
        for (int i = 0; i < r; ++i) {
            const Precision xi = x[i];
            if (xi != Precision(0)) {
                col_entries.emplace_back(i, xi);
            }
        }

        // Store result (thread-safe: distinct j indices)
        out.X_cols[j] = std::move(col_entries);
    }

    return out;
}

}} // namespace fem::constraint
