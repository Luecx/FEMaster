#pragma once
/**
 * @file cb_substitute.h
 * @brief Apply single-NNZ fixes into (C,d), drop those rows, and zero/remove fixed columns.
 *
 * Steps
 * -----
 * 1) Build maps:
 *      is_fixed_col[j] = 1 if column j is fixed by a single-NNZ row.
 *      fixed_val[j]    = d_i / a_ij for that row (or 0 if homogeneous).
 *      keep_row[i]     = 0 for the single-NNZ rows (we’ll drop them for QR).
 * 2) Update d ← d − Σ_j C(:,j) * fixed_val[j] (sparse column walks).
 * 3) Rebuild C without fixed columns (entire columns skipped).
 *
 * Notes
 * -----
 * - The output C has the same shape m×n (with fixed columns empty), which is fine
 *   because later the compression stage will drop all-zero columns and build C_use.
 */

#include <Eigen/Sparse>
#include <vector>
#include "cb_types.h"

namespace fem { namespace constraint {

struct SubstituteResult {
    SparseMatrix   C;                ///< Updated C with fixed columns zeroed (to be compressed later)
    DynamicVector  d;                ///< Updated right-hand side
    std::vector<char> keep_row;      ///< Rows to keep for QR (0 for single-NNZ rows)
    std::vector<char> is_fixed_col;  ///< Per-column flags
    std::vector<Precision> fixed_val;///< Per-column fixed values
};

/**
 * @brief Apply fixed-DOF substitutions, drop those rows, and rebuild C without fixed columns.
 *
 * @param C_in    Original constraint matrix C.
 * @param d_in    Original RHS vector d (may be empty).
 * @param singles List of single-NNZ rows from scan_single_nnz_rows().
 * @return SubstituteResult with updated (C,d) and bookkeeping vectors.
 */
SubstituteResult substitute_fixed(const SparseMatrix& C_in,
                                  const DynamicVector& d_in,
                                  const std::vector<SimpleRow>& singles);

}} // namespace fem::constraint
