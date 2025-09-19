/******************************************************************************
* @file preprocess.h
 * @brief Preprocess constraints: detect direct DOF fixes, apply substitutions,
 *        and compress zero columns with row filtering.
 *
 * Steps:
 *  1) Scan single-NNZ rows   → direct fixed DOFs u[j] = d/a.
 *  2) Substitute into (C,d)   → d := d − C(:,j)*ufix; zero fixed columns.
 *  3) Keep rows not handled by (1); compress zero columns → C_use, used list.
 *
 * All operations are sparse; no dense QR or dense intermediates.
 *
 * @date 19.09.2025
 * @author Finn
 ******************************************************************************/

#pragma once
#include "../../core/types_eig.h"
#include "../constraint_set.h"
#include <vector>

namespace fem { namespace constraint {

struct PreprocessInput {
    SparseMatrix   C;       ///< original C (will not be modified)
    DynamicVector  d;       ///< original d (empty or size m)
    int            n;       ///< number of DOFs (columns)
    int            m;       ///< number of equations (rows)
    bool           homogeneous; ///< true if d == 0
};

struct PreprocessOutput {
    // Kept rows mask (for rank/redundancy stats)
    std::vector<char> keep_row;

    // Column filtering result
    std::vector<int>  used;       ///< old column indices kept in C_use
    Eigen::VectorXi   old2new;    ///< size n, -1 if dropped

    // Fixed columns from single-NNZ rows
    std::vector<char>     is_fixed_col; ///< size n
    std::vector<Precision> fixed_val;   ///< size n

    // Modified system for QR
    SparseMatrix  C_use;  ///< m' x n_use (rows are original indices; cols compacted)
    DynamicVector d_mod;  ///< size m (matches original row indexing)

    // Counts for quick reporting
    int n_use = 0;        ///< number of used columns (nonzero in kept rows)
    int n_fixed = 0;      ///< number of directly fixed columns
};

PreprocessOutput preprocess_constraints(const PreprocessInput& in);

}} // namespace
