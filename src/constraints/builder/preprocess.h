/******************************************************************************
 * @file preprocess.h
 * @brief Declares preprocessing steps for constraint assembly.
 *
 * Preprocessing detects directly fixable DOFs, substitutes them into the
 * constraint system, and compacts the remaining rows and columns.
 *
 * @see src/constraints/builder/preprocess.cpp
 * @see src/constraints/constraint_set.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../../core/types_eig.h"
#include "../constraint_set.h"

#include <vector>

namespace fem {
namespace constraint {

/******************************************************************************
 * @struct PreprocessInput
 * @brief Raw constraint data supplied to the preprocessing stage.
 ******************************************************************************/
struct PreprocessInput {
    SparseMatrix C;        ///< Original constraint matrix.
    DynamicVector d;       ///< Original right-hand side (may be empty).
    int n = 0;             ///< Number of DOFs (columns).
    int m = 0;             ///< Number of equations (rows).
    bool homogeneous = true; ///< True if `d` is numerically zero.
};

/******************************************************************************
 * @struct PreprocessOutput
 * @brief Result of preprocessing including compacted matrices and metadata.
 ******************************************************************************/
struct PreprocessOutput {
    std::vector<char> keep_row; ///< Mask indicating rows retained after filtering.

    std::vector<int> used;      ///< Mapping of retained columns to original indices.
    Eigen::VectorXi old2new;    ///< Maps original columns to compacted indices (`-1` if dropped).

    std::vector<char> is_fixed_col;    ///< Flags indicating directly fixed DOFs.
    std::vector<Precision> fixed_val;  ///< Values assigned to fixed DOFs.

    SparseMatrix C_use; ///< Compacted constraint matrix.
    DynamicVector d_mod;///< Modified right-hand side aligned to original rows.

    int n_use = 0;   ///< Number of columns kept for QR.
    int n_fixed = 0; ///< Number of columns fixed directly.
};

/******************************************************************************
 * @brief Preprocesses the constraint system before QR factorisation.
 *
 * @param input Input constraint data.
 * @return PreprocessOutput Filtered matrices and metadata.
 ******************************************************************************/
PreprocessOutput preprocess_constraints(const PreprocessInput& input);

} // namespace constraint
} // namespace fem
