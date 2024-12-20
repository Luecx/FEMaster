/******************************************************************************
* @file reduce_mat.h
 * @brief Provides a function that reduces a sparse matrix by removing rows
 * that correspond to non-`NaN` values in a vector.
 *
 * The `reduce_mat` function reduces a sparse matrix by extracting rows
 * corresponding to `NaN` values in the given vector and removing rows that
 * correspond to non-`NaN` values.
 *
 * @param matrix The input sparse matrix.
 * @param b The input vector that controls which rows of the matrix to discard.
 * @return SparseMatrix The reduced sparse matrix, containing only rows
 * corresponding to `NaN` values in `b`.
 *
 * @date Created on 28.08.2024
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>

namespace fem { namespace mattools {

/******************************************************************************
 * @brief Reduces a sparse matrix by removing rows that correspond to non-`NaN`
 * values in the vector `b`.
 *
 * For each element in `b`, if the value is `NaN`, the corresponding row from
 * the sparse matrix is kept. Otherwise, it is discarded.
 *
 * @param matrix The input sparse matrix.
 * @param b The input vector that indicates which rows of the matrix to discard
 * (non-`NaN` values).
 * @return SparseMatrix The reduced sparse matrix containing only rows corresponding
 * to `NaN` values in `b`.
 ******************************************************************************/
SparseMatrix reduce_mat_to_mat(const SparseMatrix& matrix, const DynamicVector& b);

} } // namespace fem::mattools
