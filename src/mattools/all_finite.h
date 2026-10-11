/**
 * @file all_finite.h
 * @brief Declares a finiteness check for sparse system matrices.
 *
 * The matrix tools layer inspects numerical coefficients independently of their
 * physical meaning. Callers remain responsible for reporting invalid matrices
 * and deciding whether an assembly or solution step may proceed. Sparse storage
 * and scalar precision are defined by the core Eigen type aliases.
 *
 * @see all_finite
 * @author Finn Eggers
 * @date 11.10.2026
 */

#pragma once

#include "../core/types_eig.h"

namespace fem {
namespace mattools {

// Check every stored coefficient for NaN and positive/negative infinity.
// Implicit zeros are finite; empty matrices therefore return true. Support
// rectangular matrices and compressed/uncompressed storage without modification.
// Return false at the first nonfinite entry; leave error reporting to the caller.
bool all_finite(const SparseMatrix& matrix);

} // namespace mattools
} // namespace fem
