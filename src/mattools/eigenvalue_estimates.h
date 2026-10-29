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

#include "../constraints/transformer/constraint_map.h"
#include "../core/types_eig.h"
#include "../data/field.h"

namespace fem {
namespace mattools {

// Return sum(A_ii) / sum(B_ii) for positive semidefinite modal operators.
// This characteristic scale need not approximate the lowest eigenvalue.
// Return NaN for incompatible dimensions or a nonfinite/nonpositive scale.
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B);

// Estimate the smallest positive lambda in A q = lambda B q using at most
// 18 quadratic displacement fields with compatible beam/shell rotations.
// A must be positive on the flexible trial space. B may be indefinite.
// Return NaN when no finite positive Ritz value can be determined.
Precision estimate_eigenvalue_ritz(const SparseMatrix&              A,
                                   const SparseMatrix&              B,
                                   const model::Field&              positions,
                                   const IndexMatrix&               active_dof_idx_mat,
                                   const constraint::ConstraintMap& map);

} // namespace mattools
} // namespace fem
