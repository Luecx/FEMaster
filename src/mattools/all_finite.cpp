/**
 * @file all_finite.cpp
 * @brief Implements finiteness checks over sparse matrix coefficients.
 *
 * The matrix tools layer traverses stored coefficients without materializing
 * implicit zeros or changing the sparse storage state. Error messages and
 * solver-specific validation remain responsibilities of the calling subsystem.
 *
 * @see all_finite
 * @author Finn Eggers
 * @date 11.10.2026
 */

#include "all_finite.h"

#include <cmath>

namespace fem {
namespace mattools {

/**
 * Determines whether every sparse matrix coefficient is finite.
 *
 * Only stored entries need inspection because implicit zeros are finite. The
 * traversal stops at the first NaN or positive/negative infinity. Explicitly
 * stored zeros are inspected as well. Both compressed and uncompressed Eigen
 * storage are supported, and the input matrix remains unchanged.
 *
 * @param matrix Sparse matrix of arbitrary, including rectangular, shape.
 * @return True if all coefficients are finite; also true for an empty matrix.
 */
bool all_finite(const SparseMatrix& matrix) {
    // Traverse stored column entries; implicit zeros require no validation.
    for (Eigen::Index col = 0; col < matrix.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator entry(matrix, col); entry; ++entry) {
            if (!std::isfinite(entry.value())) {
                return false;
            }
        }
    }

    return true;
}

} // namespace mattools
} // namespace fem
