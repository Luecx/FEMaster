/******************************************************************************
* @file lump_matrix.h
 * @brief Provides a function that computes a lumped mass matrix by summing
 * each row of a sparse matrix and placing the result on the diagonal.
 *
 * The `lump_matrix` function computes a lumped mass matrix by taking a sparse
 * matrix as input, summing the values in each row, and creating a diagonal
 * matrix where the diagonal elements are the row sums.
 *
 * @param matrix The input sparse matrix.
 * @return SparseMatrix The lumped mass matrix as a diagonal sparse matrix.
 *
 * @date Created on 28.08.2024
 ******************************************************************************/

#pragma once

#include "../core/types.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace fem { namespace mattools {

/******************************************************************************
 * @brief Computes a lumped mass matrix by summing each row of the input sparse
 * matrix and placing the result on the diagonal.
 *
 * The lumped mass matrix is constructed as a diagonal matrix where the diagonal
 * entries are the sum of each row of the input sparse matrix.
 *
 * @param matrix The input sparse matrix.
 * @return SparseMatrix The lumped mass matrix as a diagonal sparse matrix.
 ******************************************************************************/
SparseMatrix lump_matrix(const SparseMatrix& matrix);

} } // namespace fem::mattools
