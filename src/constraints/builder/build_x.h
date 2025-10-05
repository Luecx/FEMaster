/******************************************************************************
 * @file build_x.h
 * @brief Declares helpers to extract sparse `X` columns from QR factors.
 *
 * The routine computes the blocks of `X = -R11^{-1} R12` by performing
 * upper-triangular back substitutions on the sparse QR factor.
 *
 * @see src/constraints/builder/build_x.cpp
 * @see src/constraints/builder/factorize_qr.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../../core/types_eig.h"

#include <utility>
#include <vector>

namespace fem {
namespace constraint {

/******************************************************************************
 * @struct XCols
 * @brief Stores sparse column data for the `X` matrix.
 ******************************************************************************/
struct XCols {
    std::vector<std::vector<std::pair<int, Precision>>> cols; ///< Column-wise entries of `X`.
};

/******************************************************************************
 * @brief Builds the columns of `X = -R11^{-1} R12` using back substitution.
 *
 * @param R Upper-triangular factor from sparse QR.
 * @param r Numerical rank of the constraint matrix.
 * @return XCols Sparse column representation of `X`.
 ******************************************************************************/
XCols build_X_cols_from_R(const SparseMatrix& R, int r);

} // namespace constraint
} // namespace fem
