/**
 * @file max_abs_diagonal.h
 * @brief Declares the maximum absolute diagonal coefficient of a sparse matrix.
 *
 * The matrix tools layer provides a numerical diagonal magnitude independent
 * of stiffness, mass or constraint semantics. Callers choose any fallback scale
 * and regularization based on this magnitude. Matrix storage and scalar
 * precision are supplied by the core Eigen type aliases.
 *
 * @see max_abs_diagonal
 * @author Finn Eggers
 * @date 11.10.2026
 */

#pragma once

#include "../core/types_eig.h"

namespace fem {
namespace mattools {

// Return max_i |matrix(i,i)|, treating missing diagonal entries as zero.
// An empty/all-zero diagonal returns zero; rectangular matrices use their main
// diagonal. Ignore off-diagonal values and leave the source matrix unchanged.
// Return NaN for any diagonal NaN, otherwise positive infinity for a diagonal
// infinity. The caller chooses validation, fallback scales and regularization.
Precision max_abs_diagonal(const SparseMatrix& matrix);

} // namespace mattools
} // namespace fem
