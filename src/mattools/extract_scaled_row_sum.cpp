/**
* @file extract_scaled_row_sum.cpp
 * @brief Implements a function that extracts rows of a sparse matrix, scales
 * them by corresponding entries in a vector, and sums the result. `NaN` values
 * in the vector are ignored.
 *
 * The function multiplies each row of the sparse matrix by the non-`NaN` values
 * in the vector, then sums the scaled rows into a final output vector.
 *
 * The sparse matrix representation improves performance for larger data sets.
 *
 * @date Created on 28.08.2024
 */

#include "extract_scaled_row_sum.h"

namespace fem { namespace mattools {

DynamicVector extract_scaled_row_sum(const SparseMatrix& matrix, const DynamicMatrix& scale_vector) {
    // ensure square matrix
    logging::error(matrix.rows() == matrix.cols(), "Matrix must be square.");
    logging::error(matrix.rows() == scale_vector.size(), "Matrix row count must match the size of the scale_vector.");

    // Initialize the result vector with the appropriate size (matching the number of columns in the matrix)
    DynamicVector result = DynamicVector::Zero(matrix.cols(), 1);

    // Iterate through the scale_vector and corresponding rows in the sparse matrix
    for (int i = 0; i < scale_vector.size(); ++i) {
        // Check if the current value in the scale_vector is not NaN
        if (!std::isnan(scale_vector(i))) {

            if (scale_vector(i) == 0) {
                continue;
            }

            // Extract the corresponding row from the sparse matrix
            SparseMatrix::InnerIterator it(matrix, i);

            // Iterate over the non-zero elements in the sparse row
            while (it) {
                // Add the scaled value of the row entry to the result
                result(it.row()) += it.value() * scale_vector(i);
                ++it;
            }
        }
    }

    return result;
}

} } // namespace fem::mattools
