/******************************************************************************
 * @file extract_scaled_row_sum.h
 * @brief Provides a function that extracts rows of a sparse matrix, scales
 * them by corresponding entries in a vector, and sums the result. The vector
 * may contain `NaN` values, which are ignored.
 *
 * Given a sparse matrix and a vector, this function multiplies each row of
 * the matrix by the corresponding non-`NaN` value in the vector and sums
 * the results to produce a single output vector.
 *
 * The matrix uses a sparse representation to handle larger data efficiently.
 *
 * @param matrix The input sparse matrix with rows to be scaled and summed.
 * @param scale_vector A vector of scalars, some of which may be `NaN`.
 * Non-`NaN` entries are used to scale the corresponding rows of the matrix.
 * @return DynamicMatrix The result of summing all scaled rows into a single
 * output vector.
 *
 * @date Created on 28.08.2024
 ******************************************************************************/

#pragma once

#include "../core/types.h"
#include "../core/logging.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>  // For std::isnan

namespace fem { namespace mattools {

/******************************************************************************
 * @brief Extracts non-`NaN` entries from a vector, multiplies corresponding
 * rows of a sparse matrix by the vector's values, and sums the resulting scaled rows.
 *
 * Each non-`NaN` entry in the vector is used to scale the corresponding row
 * in the sparse matrix. The scaled rows are then summed into a single output vector.
 *
 * The input matrix uses a sparse representation to handle large data more efficiently.
 *
 * @param matrix The input sparse matrix with rows to be scaled and summed.
 * @param scale_vector A vector of scalars, where some values may be `NaN`.
 * @return DynamicMatrix The resulting vector of the sum of all scaled rows.
 ******************************************************************************/
DynamicVector extract_scaled_row_sum(const SparseMatrix& matrix, const DynamicMatrix& scale_vector);

} } // namespace fem::mattools
