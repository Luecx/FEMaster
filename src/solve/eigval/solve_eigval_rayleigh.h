#pragma once

#include "../../constraints/transformer/constraint_map.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"

namespace fem::solver {

/**
 * @brief Estimate a low flexible eigenvalue scale from a quadratic Rayleigh-Ritz subspace.
 *
 * The estimator uses the six scalar monomials x², y², z², xy, xz and yz on
 * normalized nodal coordinates and applies each in the three translational
 * directions. Compatible nodal rotations are added from 1/2 curl(u) for
 * beam/shell rotational DOFs. Admissible rigid-body modes are detected from
 * the reduced stiffness matrix, removed in the B inner product, and the
 * remaining trial vectors are B-orthonormalized. The smallest positive Ritz
 * value of the projected generalized eigenproblem is returned.
 */
Precision estimate_lambda_scale_from_geometry(
    const SparseMatrix&              A,
    const SparseMatrix&              B,
    const model::Field&              positions,
    const IndexMatrix&               active_dof_idx_mat,
    const constraint::ConstraintMap& map
);

} // namespace fem::solver
