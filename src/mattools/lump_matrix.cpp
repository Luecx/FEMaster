/******************************************************************************
* @file lump_matrix.cpp
 * @brief Implementation of the function that computes a lumped mass matrix by
 * summing each row of a sparse matrix and placing the result on the diagonal.
 *
 * @date Created on 28.08.2024
 ******************************************************************************/

#include "lump_matrix.h"

namespace fem { namespace mattools {

SparseMatrix lump_matrix(const SparseMatrix& matrix) {
    // Number of rows in the input matrix
    int rows = matrix.rows();

    // Create a diagonal sparse matrix to store the lumped mass matrix
    SparseMatrix lumped_matrix(rows, rows);

    // Vector to store the row sums
    DynamicVector row_sums(rows);

    // Iterate over each row and sum the non-zero elements
    for (int i = 0; i < rows; ++i) {
        row_sums[i] = matrix.row(i).sum();
    }

    // Fill the diagonal of the sparse matrix with the row sums
    lumped_matrix.reserve(rows);
    for (int i = 0; i < rows; ++i) {
        if (row_sums[i] != 0) {
            lumped_matrix.insert(i, i) = row_sums[i];
        }
    }

    // Compress the sparse matrix to optimize storage
    lumped_matrix.makeCompressed();

    return lumped_matrix;
}

} } // namespace fem::mattools
