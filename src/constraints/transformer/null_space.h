/**
 * @file null_space.h
 * @brief Declares construction of affine null-space constraint maps.
 *
 * @see src/constraints/transformer/null_space.cpp
 * @see src/constraints/transformer/constraint_map.h
 * @author Finn Eggers
 * @date 14.07.2026
 */

#pragma once

#include "constraint_build_report.h"
#include "constraint_map.h"
#include "constraint_system.h"

#include <utility>

namespace fem {
namespace constraint {

/**
 * @brief Controls numerical rank detection and feasibility checks.
 */
struct NullSpaceOptions {
    Precision rank_tolerance        = Precision(1e-12); ///< Relative threshold applied to the diagonal of `R`.
    Precision feasibility_tolerance = Precision(1e-10); ///< Relative residual tolerance scaled by `max(||d||, 1)`.
};

/**
 * @brief Builds the affine transformation `u = u_p + T q` from `C u = d`.
 *
 * Directly prescribed DOFs are substituted before a column-pivoted sparse QR
 * factorization determines the remaining dependent and independent DOFs.
 *
 * @param system Assembled constraint system over active full-system DOFs.
 * @param options Numerical tolerances for rank and feasibility detection.
 * @return Null-space map and diagnostics for the assembled constraint system.
 */
std::pair<ConstraintMap, ConstraintBuildReport>
build_null_space(const ConstraintSystem& system, const NullSpaceOptions& options = {});

} // namespace constraint
} // namespace fem
