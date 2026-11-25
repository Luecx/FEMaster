/**
 * @file extract_connected_rows.h
 * @brief Declares extraction of rows coupled to seed rows through shared columns.
 *
 * Two rows are coupled when they share a nonzero column. Starting from a set of
 * seed rows, the matrix tools layer collects the transitive closure of this
 * coupling and returns the corresponding row block. Callers use the returned
 * row indices to map local results back to the original numbering.
 *
 * @see extract_connected_rows
 * @author Finn Eggers
 * @date 25.11.2026
 */

#pragma once

#include "../core/types_eig.h"

#include <vector>

namespace fem {
namespace mattools {

// Return the rows of matrix reachable from the seed rows through shared nonzero
// columns. The result keeps all columns, renumbers the selected rows locally in
// ascending original order and writes their original indices to rows.
// Duplicate seeds are allowed; an empty seed list yields an empty matrix.
SparseMatrix extract_connected_rows(const SparseMatrix&       matrix,
                                    const std::vector<Index>& seeds,
                                    std::vector<Index>&       rows);

} // namespace mattools
} // namespace fem
