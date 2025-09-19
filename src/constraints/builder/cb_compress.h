#pragma once
/**
 * @file cb_compress.h
 * @brief Column compression of a sparse constraint matrix with optional row filtering.
 *
 * Purpose
 * -------
 * Build a column-compressed submatrix `C_use` by:
 *   1) (optionally) dropping rows where `keep_row[i] == 0`
 *   2) removing all columns that are entirely zero on the kept rows
 *   3) recording the surviving column indices and an index remapping
 *
 * This is used before SparseQR so we only factor columns that actually
 * participate in the remaining (kept) equations.
 *
 * Interface
 * ---------
 *   compress_zero_columns_with_row_filter(C, keep_row, used_cols, old2new, C_use)
 *
 * Inputs:
 *   - C         : m×n sparse matrix in column-major compressed storage.
 *   - keep_row  : size m. If empty, all rows are kept.
 *                 Otherwise, rows with keep_row[i] == 0 are ignored (treated as removed).
 *
 * Outputs:
 *   - used_cols : list of global column indices (in original [0..n-1]) that have
 *                 at least one nonzero in the kept rows. The order follows the
 *                 original column order.
 *   - old2new   : length n; old2new[j] = new column index in C_use if column j survived, or -1 otherwise.
 *   - C_use     : m×k sparse matrix (k = used_cols.size()), formed by copying the
 *                 surviving columns of C and discarding entries in dropped rows.
 *                 Row indices are NOT renumbered; C_use has the same number of rows (m).
 *
 * Guarantees & Notes
 * ------------------
 * - Complexity: O(nnz(C)) to scan + copy (with cheap checks). Memory is linear in nnz(C_use).
 * - Stability: we do exact-zero tests only; no numeric thresholding here.
 * - If `keep_row` is empty, this reduces to a plain “drop all-zero columns” compression.
 * - `C_use` is returned compressed (makeCompressed()).
 * - `used_cols.size()` equals the number of non-empty columns on the kept rows.
 * - `old2new[j] == -1` iff column j becomes empty after filtering the kept rows.
 */

#include <Eigen/Sparse>
#include <vector>
#include "../../core/types_eig.h"

namespace fem { namespace constraint {

/**
 * @brief Build `C_use` by removing columns of `C` that are zero on the kept rows.
 *
 * @param C         Input m×n sparse matrix.
 * @param keep_row  Row-keep mask of length m. If empty, all rows are kept.
 * @param used_cols Output: global column indices that survive (size k).
 * @param old2new   Output: length n; maps old column j -> [0..k-1] or -1 if dropped.
 * @param C_use     Output: m×k sparse matrix consisting of surviving columns of C,
 *                           with entries in dropped rows omitted.
 */
void compress_zero_columns_with_row_filter(const SparseMatrix& C,
                                           const std::vector<char>& keep_row,
                                           std::vector<int>&   used_cols,
                                           Eigen::VectorXi&    old2new,
                                           SparseMatrix&       C_use);

}} // namespace fem::constraint
