/**
 * @file append_sparse_block.cpp
 * @brief Implements offset and scaling of sparse coefficients into triplets.
 *
 * The matrix tools layer appends entries according to
 * (i,j,A_ij) -> (row_offset+i,col_offset+j,factor*A_ij).
 * It validates representability in the core Triplet storage index type before
 * insertion. Destination dimensions, physical block interpretation and duplicate
 * accumulation remain responsibilities of the caller and Eigen construction.
 *
 * @see append_sparse_block
 * @author Finn Eggers
 * @date 11.10.2026
 */

#include "append_sparse_block.h"

#include "../core/logging.h"

#include <limits>

namespace fem {
namespace mattools {

/**
 * Appends a scaled sparse matrix block to an existing triplet buffer.
 *
 * Each stored source coefficient A(i,j) contributes factor*A(i,j) at
 * (row_offset+i,col_offset+j). Independent offsets support rectangular blocks.
 * Existing buffer entries remain intact, and duplicate coordinates are retained
 * for summation by a subsequent Eigen setFromTriplets call.
 *
 * Offset and block extent validation precedes insertion to prevent overflow
 * when converting FEMaster's unsigned Index to Eigen's Triplet::StorageIndex.
 * The destination matrix is not supplied, so its dimensions must be checked by
 * the caller. The caller also owns buffer capacity planning across all blocks.
 *
 * Explicit zeros and nonfinite source coefficients are preserved under scalar
 * multiplication. No pruning or finiteness checks are performed: multiplication
 * by zero follows ordinary floating-point arithmetic, including for NaN/Inf.
 * An empty block appends nothing. The source remains unchanged and may use
 * compressed or uncompressed storage. Buffer allocation may throw; allocation
 * failure during insertion can leave a partially appended block.
 *
 * @param triplets Caller-owned buffer receiving additional contributions.
 * @param matrix Source block with its original row and column numbering.
 * @param row_offset Nonnegative destination offset for source rows.
 * @param col_offset Nonnegative destination offset for source columns.
 * @param factor Scalar multiplier, including negative or zero values.
 */
void append_sparse_block(TripletList&        triplets,
                         const SparseMatrix& matrix,
                         Index               row_offset,
                         Index               col_offset,
                         Precision           factor) {
    // Validate offsets before subtraction, then verify the largest block indices.
    using StorageIndex = Triplet::StorageIndex;
    const Index maximum_index = static_cast<Index>(std::numeric_limits<StorageIndex>::max());
    logging::error(row_offset <= maximum_index && col_offset <= maximum_index,
        "append_sparse_block: offsets exceed triplet storage index range");
    logging::error(matrix.rows() == 0 || static_cast<Index>(matrix.rows() - 1) <= maximum_index - row_offset,
        "append_sparse_block: row extent exceeds triplet storage index range");
    logging::error(matrix.cols() == 0 || static_cast<Index>(matrix.cols() - 1) <= maximum_index - col_offset,
        "append_sparse_block: column extent exceeds triplet storage index range");

    // Embed stored local entries into the destination coordinates. Leave
    // duplicate summation and sparse matrix construction to the caller.
    for (Eigen::Index col = 0; col < matrix.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator entry(matrix, col); entry; ++entry) {
            const StorageIndex row    = static_cast<StorageIndex>(row_offset + static_cast<Index>(entry.row()));
            const StorageIndex column = static_cast<StorageIndex>(col_offset + static_cast<Index>(entry.col()));
            triplets.emplace_back(row, column, factor * entry.value());
        }
    }
}

} // namespace mattools
} // namespace fem
