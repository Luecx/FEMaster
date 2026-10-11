/**
 * @file append_sparse_block.h
 * @brief Declares insertion of scaled sparse blocks into a triplet buffer.
 *
 * The matrix tools layer maps local matrix entries into a larger algebraic
 * matrix through row and column offsets. Callers retain responsibility for
 * destination dimensions, buffer reservation, block layout and final sparse
 * construction. Eigen combines duplicate contributions during setFromTriplets.
 *
 * @see append_sparse_block
 * @author Finn Eggers
 * @date 11.10.2026
 */

#pragma once

#include "../core/types_eig.h"

namespace fem {
namespace mattools {

// Append each stored matrix(i,j) as factor * matrix(i,j) at destination position
// (row_offset+i, col_offset+j). Support rectangular blocks and preserve existing
// triplets, duplicate coordinates and explicitly stored zeros. Empty blocks add
// nothing; nonfinite coefficients follow ordinary scalar multiplication.
// Validate that offsets and block indices fit SparseMatrix::StorageIndex before
// insertion. The caller reserves buffer capacity, checks destination dimensions
// and constructs the final sparse matrix, where Eigen can sum duplicate entries.
void append_sparse_block(TripletList&        triplets,
                         const SparseMatrix& matrix,
                         Index               row_offset = 0,
                         Index               col_offset = 0,
                         Precision           factor     = Precision(1));

} // namespace mattools
} // namespace fem
