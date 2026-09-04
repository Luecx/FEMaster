#pragma once

#include "../../constraints/transformer/constraint_map.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"

namespace fem::solver {

/**
 * @brief Estimate a low flexible eigenvalue scale from quadratic geometry-based trial fields.
 *
 * The estimator uses the six scalar monomials x², y², z², xy, xz and yz on
 * normalized nodal coordinates and applies each independently in the three
 * translational directions. Admissible rigid-body translations and rotations
 * are detected from the reduced stiffness matrix and removed from every trial
 * vector in the B inner product before evaluating the Rayleigh quotient. The
 * smallest positive quotient is returned.
 */
Precision estimate_lambda_scale_from_geometry(
    const SparseMatrix&              A,
    const SparseMatrix&              B,
    const model::Field&              positions,
    const IndexMatrix&               active_dof_idx_mat,
    const constraint::ConstraintMap& map
);

} // namespace fem::solver
