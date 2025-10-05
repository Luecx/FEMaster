/**
 * @file invariants.h
 * @brief Declares debug-time invariant checks for constraint maps.
 *
 * The helper verifies core identities such as `C T = 0` to catch programming
 * errors during development.
 *
 * @see src/constraints/builder/invariants.cpp
 * @see src/constraints/constraint_map.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../../core/types_eig.h"
#include "../constraint_map.h"

namespace fem {
namespace constraint {

/**
 * @brief Validates basic invariants of the assembled constraint map.
 *
 * @param map Constraint map to test.
 * @param C Original constraint matrix.
 */
void check_invariants(const ConstraintMap& map, const SparseMatrix& C);

} // namespace constraint
} // namespace fem
