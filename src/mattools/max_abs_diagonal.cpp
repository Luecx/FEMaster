/**
 * @file max_abs_diagonal.cpp
 * @brief Implements sparse diagonal magnitude evaluation.
 *
 * The matrix tools layer computes max_i |A(i,i)| without allocating a dense
 * diagonal vector. It preserves nonfinite diagonal values for caller-side
 * validation. Physical scaling choices and regularization remain outside this
 * numerical operation.
 *
 * @see max_abs_diagonal
 * @author Finn Eggers
 * @date 11.10.2026
 */

#include "max_abs_diagonal.h"

#include <algorithm>
#include <cmath>

namespace fem {
namespace mattools {

/**
 * Computes the maximum absolute value on the main diagonal.
 *
 * Stored diagonal coefficients are accumulated column by column. Missing
 * coefficients are implicit zeros, so initializing the maximum to zero also
 * handles an empty diagonal. For rectangular matrices, only entries with
 * matching row and column indices contribute. Off-diagonal values, including
 * nonfinite ones, do not affect this diagonal-only measurement.
 *
 * A diagonal NaN is returned explicitly rather than silently ignored by a
 * maximum comparison. A diagonal infinity produces positive infinity unless a
 * diagonal NaN is encountered. Both sparse storage states are supported, and
 * no input coefficients or metadata are modified.
 *
 * @param matrix Sparse matrix whose main diagonal is measured.
 * @return Maximum diagonal magnitude, zero for an empty/all-zero diagonal,
 *         NaN for a diagonal NaN, or positive infinity for a diagonal infinity.
 */
Precision max_abs_diagonal(const SparseMatrix& matrix) {
    // Start from the magnitude of implicit zeros and inspect stored diagonals.
    Precision maximum = Precision(0);
    for (Eigen::Index col = 0; col < matrix.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator entry(matrix, col); entry; ++entry) {
            if (entry.row() != entry.col()) {
                continue;
            }

            // Propagate NaN explicitly; std::max would otherwise conceal it.
            const Precision magnitude = std::abs(entry.value());
            if (std::isnan(magnitude)) {
                return magnitude;
            }
            maximum = std::max(maximum, magnitude);
        }
    }

    return maximum;
}

} // namespace mattools
} // namespace fem
