/**
 * @file regularise_stiffness.cpp
 * @brief Implements row-based regularisation of sparse stiffness matrices.
 *
 * The implementation identifies numerically weak rows through their absolute
 * coefficient sums and adds positive diagonal contributions where required.
 *
 * A representative stiffness scale is first obtained from the mean absolute
 * value of the non-zero diagonal coefficients. If no non-zero diagonal entry
 * exists, the mean absolute value of all non-zero matrix coefficients is used
 * as a fallback.
 *
 * Each row whose absolute coefficient sum lies below the resulting
 * matrix-relative threshold receives an additional diagonal stiffness. This
 * prevents empty or weakly constrained degrees of freedom from producing
 * exactly singular rows during the subsequent linear solve.
 *
 * The routine performs only numerical matrix stabilisation. It does not alter
 * the active degree-of-freedom set, constraint equations or load vectors.
 *
 * @see regularise_stiffness
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#include "regularise_stiffness.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {

/**
 * Adds positive diagonal regularisation to numerically weak matrix rows.
 *
 * For every row, the routine computes the absolute coefficient sum
 *
 *     row_sum(i) = sum_j |K(i, j)|.
 *
 * A representative global stiffness scale is determined from the mean
 * absolute value of all non-zero diagonal entries. If the matrix does not
 * contain a non-zero diagonal entry, the mean absolute value of all non-zero
 * matrix coefficients is used instead.
 *
 * The minimum admissible row stiffness is defined as
 *
 *     max(alpha * stiffness_scale, machine_epsilon).
 *
 * For every row below this threshold, the missing stiffness is added to the
 * corresponding diagonal entry. All additions are collected as sparse
 * triplets and applied to the matrix in one operation.
 *
 * No modification is performed when:
 *
 * - `alpha` is non-positive,
 * - the matrix contains no rows,
 * - no finite positive stiffness scale can be determined,
 * - every row already satisfies the minimum stiffness threshold.
 *
 * @param matrix Sparse stiffness matrix modified directly by the routine.
 * @param alpha Relative regularisation factor applied to the global stiffness
 *              scale. A non-positive value disables regularisation.
 *
 * @return Number of matrix rows that received an additional diagonal term.
 */
int regularise_stiffness(SparseMatrix& matrix, Precision alpha) {
    // Treat a non-positive factor as disabled and avoid processing an empty
    // stiffness matrix
    if (alpha <= Precision(0) || matrix.rows() == 0) {
        return 0;
    }

    // Accumulate the data required for the global stiffness scale and the
    // absolute coefficient sum of every matrix row
    Precision     diagonal_sum   = Precision(0);
    Index         diagonal_count = 0;
    Precision     entry_sum      = Precision(0);
    Index         entry_count    = 0;
    DynamicVector row_sum        = DynamicVector::Zero(matrix.rows());

    // Traverse the sparse matrix once. Eigen stores the matrix by outer
    // indices, which are columns for the default column-major SparseMatrix.
    for (Index col = 0; col < static_cast<Index>(matrix.outerSize()); ++col) {
        for (SparseMatrix::InnerIterator it(matrix, col); it; ++it) {
            const Precision absolute_value = std::abs(it.value());

            // Measure the total absolute stiffness associated with the row
            row_sum(it.row()) += absolute_value;

            // Collect all non-zero entries as a fallback scale in case the
            // matrix does not contain usable diagonal coefficients
            if (absolute_value > Precision(0)) {
                entry_sum += absolute_value;
                ++entry_count;
            }

            // Prefer diagonal coefficients for the representative stiffness
            // scale because they most directly describe individual DOF
            // stiffnesses
            if (it.row() == it.col() && absolute_value > Precision(0)) {
                diagonal_sum += absolute_value;
                ++diagonal_count;
            }
        }
    }

    // Determine the representative matrix stiffness. The mean absolute
    // diagonal is preferred, while the mean absolute matrix entry serves as a
    // fallback for matrices without non-zero diagonal coefficients.
    Precision stiffness_scale = Precision(0);

    if (diagonal_count > 0) {
        stiffness_scale = diagonal_sum / static_cast<Precision>(diagonal_count);
    } else if (entry_count > 0) {
        stiffness_scale = entry_sum / static_cast<Precision>(entry_count);
    }

    // Without a finite positive reference scale, no meaningful relative
    // regularisation stiffness can be constructed
    if (stiffness_scale <= Precision(0) || !std::isfinite(stiffness_scale)) {
        return 0;
    }

    // Prevent the threshold itself from becoming exactly zero through
    // underflow when alpha or the representative stiffness is very small
    const Precision min_row_stiffness = std::max(
        alpha * stiffness_scale,
        std::numeric_limits<Precision>::epsilon()
    );

    // Collect all required diagonal additions before modifying the sparse
    // matrix. This avoids repeated sparse insertions during the row scan.
    TripletList additions{};

    for (Index row = 0; row < static_cast<Index>(matrix.rows()); ++row) {
        // Rows below the threshold receive exactly the currently missing
        // stiffness as an additional positive diagonal contribution
        if (row_sum(row) < min_row_stiffness) {
            const Precision missing_stiffness = min_row_stiffness - row_sum(row);

            additions.emplace_back(row, row, missing_stiffness);
        }
    }

    // Keep the original matrix unchanged when every row is already
    // sufficiently stiff
    if (additions.empty()) {
        return 0;
    }

    // Assemble the collected diagonal terms into a temporary sparse matrix
    SparseMatrix regularisation(matrix.rows(), matrix.cols());
    regularisation.setFromTriplets(additions.begin(), additions.end());

    // Add all regularisation terms at once and restore compressed sparse
    // storage for subsequent matrix operations and factorisation
    matrix += regularisation;
    matrix.makeCompressed();

    return static_cast<int>(additions.size());
}

} // namespace fem
