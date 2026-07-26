/**
 * @file regularise_stiffness.h
 * @brief Declares sparse stiffness-matrix regularisation.
 *
 * The regularisation utility stabilises numerically weak or empty rows of an
 * assembled sparse stiffness matrix by adding positive contributions to the
 * corresponding diagonal entries.
 *
 * The required stabilisation level is determined from a representative global
 * stiffness scale and a user-supplied relative regularisation factor. The
 * matrix is modified directly, while the function reports the number of rows
 * that received an additional diagonal contribution.
 *
 * Matrix assembly, active-degree-of-freedom selection, constraint handling and
 * linear-system solution remain responsibilities of the calling subsystem.
 *
 * @see SparseMatrix
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

// Sparse stiffness-matrix regularisation. The matrix is modified in place by
// adding positive diagonal contributions to rows whose absolute coefficient
// sum lies below the matrix-relative threshold defined by alpha.
//
// A non-positive alpha disables regularisation. The return value is the number
// of matrix rows for which a diagonal contribution was added.
int regularise_stiffness(SparseMatrix& matrix, Precision alpha);

} // namespace fem