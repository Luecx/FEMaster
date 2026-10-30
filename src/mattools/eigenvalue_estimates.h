/**
 * @file eigenvalue_estimates.h
 * @brief Declares a geometry-based Ritz estimate for structural eigenproblems.
 *
 * The matrix tools layer builds a quadratic trial space in the independent
 * coordinates of a ConstraintMap. The reduced elastic stiffness A supplies
 * the positive energy metric; B may be a mass matrix or indefinite negative
 * geometric stiffness. Constraint assembly and spectral shifts belong to the
 * calling load case. A separate diagonal-sum scale provides a geometry-independent
 * fallback for positive semidefinite modal operators. Inputs remain unchanged.
 *
 * @see constraint::ConstraintMap
 * @author Finn Eggers
 * @date 29.10.2026
 */

#pragma once

#include "../constraints/transformer/constraint_transformer.h"
#include "../core/types_eig.h"
#include "../data/field.h"

namespace fem {
namespace mattools {

// Return sum(A_ii) / sum(B_ii) for positive semidefinite modal operators.
// This characteristic scale need not approximate the lowest eigenvalue.
// Zero stiffness gives zero; the caller must still choose a nonzero shift.
// Return NaN for incompatible dimensions, a nonfinite/negative scale or nonpositive mass sum.
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B);

// Estimate the smallest positive eigenvalue using rigid and quadratic trial fields.
// The transformer must have initialized its master DOF locations beforehand.
// A must be positive semidefinite; B may be indefinite.
// Return NaN if no finite positive Ritz estimate can be determined.
Precision estimate_eigenvalue_ritz(const SparseMatrix&                      A,
                                   const SparseMatrix&                      B,
                                   const model::Field&                      positions,
                                   const constraint::ConstraintTransformer& transformer);

} // namespace mattools
} // namespace fem
