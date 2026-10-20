/**
 * @file merge_sparse_matrices.h
 * @brief Declares reduction of worker-local sparse assembly matrices.
 *
 * The matrix tools layer combines compressed partial matrices into one global
 * CSC matrix. The caller owns element evaluation, triplet insertion and input
 * compression. The merge owns structural counting, output allocation and
 * stable accumulation of overlapping entries.
 *
 * @see assemble_matrix_multithreaded
 */

#pragma once

#include "../core/types_eig.h"

#include <vector>

namespace fem::mattools {

// Sum one to 128 identically sized square CSC matrices with sorted row indices.
// Preserve explicitly stored zeros and accumulate equal coordinates in input
// order. Each output column is processed independently; inputs remain unchanged.
SparseMatrix merge_sparse_matrices(const std::vector<SparseMatrix>& matrices);

} // namespace fem::mattools
