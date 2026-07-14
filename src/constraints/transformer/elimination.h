/**
 * @file elimination.h
 * @brief Declares direct substitution of constraint equations.
 *
 * @see src/constraints/transformer/elimination.cpp
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
 * @brief Builds an affine constraint map by eliminating one dependent DOF per independent equation.
 *
 * Previously eliminated DOFs are substituted before selecting each new pivot,
 * so the resulting dependencies remain acyclic and reference only final master DOFs.
 *
 * @param system Assembled constraint system over active full-system DOFs.
 * @return Elimination map and diagnostics for the assembled constraint system.
 */
std::pair<ConstraintMap, ConstraintBuildReport> build_elimination(const ConstraintSystem& system);

} // namespace constraint
} // namespace fem
