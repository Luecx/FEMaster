/**
 * @file eigenvalue_estimates.h
 * @brief Declares numerical estimates for generalized eigenvalue problems.
 *
 * The matrix tools layer measures A x = lambda B x without assuming stiffness,
 * mass or geometric-stiffness semantics. It provides a directional Rayleigh
 * quotient, a characteristic scale from absolute diagonal sums and a lower
 * quartile of positive diagonal ratios. These quantities are not guaranteed
 * to approximate the smallest eigenvalue.
 *
 * Callers choose spectral shifts, signs, acceptance criteria and fallback
 * values. Invalid or undefined estimates return NaN; inputs remain unchanged.
 *
 * @see SparseMatrix
 * @author Finn Eggers
 * @date 27.10.2026
 */

#pragma once

#include "../core/types_eig.h"

namespace fem {
namespace mattools {

// Evaluate (x^T A x) / (x^T B x). Zero and negative quotients are valid.
// Return NaN for incompatible dimensions, nonfinite arithmetic, a zero vector
// or a denominator that vanishes through cancellation in an indefinite B.
Precision estimate_eigenvalue_rayleigh(const SparseMatrix& A, const SparseMatrix& B, const DynamicVector& x);

// Evaluate sum_i |A_ii| / sum_i |B_ii| for compatible square matrices.
// A zero numerator returns zero; null modes need not make the numerator zero.
// Return NaN for incompatible dimensions, nonfinite diagonals or a zero
// denominator. Off-diagonal coefficients do not enter this scale estimate.
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B);

// Return the lower quartile of finite positive A_ii / B_ii for A x = lambda B x.
// Retain the buckling fallback's thresholds: |A_ii| >= 1e-20 and B_ii > 1e-20.
// Skip invalid/nonpositive ratios; return NaN for incompatible square matrices
// or when no eligible diagonal entry remains. Inputs are not modified.
Precision estimate_eigenvalue_diag_ratio(const SparseMatrix& A, const SparseMatrix& B);

} // namespace mattools
} // namespace fem
