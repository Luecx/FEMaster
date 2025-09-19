#pragma once
/**
 * @file cb_build_x_sparse.h
 * @brief Sparse construction of Xd = -R11^{-1} R12 from the QR upper-triangular factor R.
 *
 * Context
 * -------
 * After a rank-revealing SparseQR of C_use, we have an upper-triangular (possibly
 * rectangular) matrix R such that (after column permutation) the leading r columns
 * form R11 and the next nm_use columns form R12. We need Xd = -R11^{-1} R12,
 * but we want to avoid any dense blocks. This step computes Xd column-by-column
 * using *sparse back-substitution* on R11 and returns it in a compact format.
 *
 * Output Format
 * -------------
 * - `X_cols`: for each QR-master j (0..nm_use-1), a sparse list of (i, value)
 *   pairs where i ∈ [0..r-1] indexes the slave rows. This is the j-th column of
 *   Xd in sparse form.
 * - `R11_diag`, `R11_rows`: a compact row-major view of the upper-triangular R11
 *   (diagonal and row-wise off-diagonal entries) that callers may reuse if needed.
 *
 * Important Notes
 * ---------------
 * - This function does *not* allocate any dense r×nm_use matrices.
 * - Each Xd column is solved independently; the implementation may parallelize
 *   across columns (e.g., with OpenMP).
 * - Only the portion of R touching rows < r is read; entries at/under the r-th
 *   row are ignored for the Xd build.
 */

#include "cb_types.h"            // RowEntry
#include <Eigen/Sparse>
#include <utility>
#include <vector>
#include "../../core/types_eig.h"

namespace fem { namespace constraint {

/// Sparse columns of Xd. X_cols[j] holds the j-th column of Xd as (row_i, value).
using XCols = std::vector<std::vector<std::pair<int, Precision>>>;

/**
 * @brief Results of the sparse Xd build.
 *
 * Members:
 * - R11_diag : diagonal of R11, length r
 * - R11_rows : for each row i in [0..r-1], a list of (k, a_ik) with k>i giving
 *              the strictly upper entries of R11(i,k). Together with R11_diag this
 *              forms a compact row-major view of R11 for back-substitution.
 * - X_cols   : sparse columns of Xd = -R11^{-1} R12; length nm_use
 */
struct XBuildOut {
    std::vector<Precision>             R11_diag;
    std::vector<std::vector<RowEntry>> R11_rows;
    XCols                               X_cols;
};

/**
 * @brief Compute Xd = -R11^{-1} R12 using sparse back-substitution on R11.
 *
 * @param R       The QR upper-triangular factor (as returned by SparseQR::matrixR()).
 *                Only the top-left r×(r+nm_use) block is consulted.
 * @param r       Numerical rank (size of R11; also number of slave rows).
 * @param nm_use  Number of QR-master columns (i.e., number of columns of R12).
 *
 * @return XBuildOut with (R11_diag, R11_rows) and sparse columns X_cols.
 *
 * Pre-conditions:
 * - 0 ≤ r ≤ min(R.rows(), R.cols()).
 * - r + nm_use ≤ R.cols().
 *
 * Complexity:
 * - Building the compact R11 view: O(nnz(R11)).
 * - Each column backsolve: proportional to the number of touched entries in that
 *   column’s dependency cone (typically small if R is sparse/triangular).
 * - Columns are independent and can be computed in parallel.
 */
XBuildOut build_X_sparse(const SparseMatrix& R, int r, int nm_use);

}} // namespace fem::constraint
