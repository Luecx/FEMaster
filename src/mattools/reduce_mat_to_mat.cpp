/**
* @file reduce_mat.cpp
 * @brief Implements a function that reduces a sparse matrix by removing rows
 * that correspond to non-`NaN` values in a vector.
 *
 * @date Created on 28.08.2024
 */

#include "reduce_mat_to_mat.h"
#include "assert.h"
#include "../core/core.h"


namespace fem { namespace mattools {

SparseMatrix reduce_mat_to_mat(const SparseMatrix& matrix, const DynamicVector& b) {
    runtime_assert(matrix.rows() == b.size(), "Matrix row count must match the size of the vector.");

    // Step 1: Create a mapping of old indices to new indices
    std::vector<int> index_map(b.size(), -1);
    int new_index = 0;

    for (int i = 0; i < b.size(); ++i) {
        if (std::isnan(b(i))) {
            index_map[i] = new_index++;
        }
    }

    // Step 2: Prepare a new reduced sparse matrix
    SparseMatrix reduced_matrix(new_index, new_index);  // Size based on number of NaNs in b
    std::vector<Eigen::Triplet<Precision>> tripletList;

    // Step 3: Iterate over the original matrix and add elements to the reduced matrix
    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            Precision value = it.value();

            // Check if both row and column should be kept
            if (index_map[row] != -1 && index_map[col] != -1) {
                tripletList.emplace_back(index_map[row], index_map[col], value);
            }
        }
    }

    // Step 4: Set the reduced matrix from the triplet list
    reduced_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    return reduced_matrix;
}

} } // namespace fem::mattools
