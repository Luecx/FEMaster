#pragma once
/**
 * @file cb_scan_simple.h
 * @brief Detect single-NNZ constraint rows (a_i u_j = d_i) to directly fix DOFs.
 *
 * For each row with exactly one nonzero at column `col` and value `a`, we can
 * deduce u[col] = d/a (or 0 if homogeneous). These DOFs are removed from QR and
 * folded into u_p directly.
 */

#include <Eigen/Sparse>
#include <vector>
#include "cb_types.h"   // SimpleRow, Precision, DynamicVector, SparseMatrix

namespace fem { namespace constraint {

/**
 * @brief Find rows of C with exactly one nonzero.
 *
 * @param C            m×n constraint matrix.
 * @param d            RHS vector (may be empty in homogeneous case).
 * @param homogeneous  If true, treat d as zero.
 * @return vector of SimpleRow{row, col, a, d_row_or_0}.
 */
std::vector<SimpleRow>
scan_single_nnz_rows(const SparseMatrix& C,
                     const DynamicVector& d,
                     bool homogeneous);

}} // namespace fem::constraint
