/**
 * @file constraint_build_report.h
 * @brief Declares diagnostic information generated while building a constraint system.
 *
 * @author Finn Eggers
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {
namespace constraint {

/**
 * @brief Summarizes the dimensions and numerical properties of a constraint system.
 */
struct ConstraintBuildReport {
    // System dimensions
    Index equations{};      ///< Number of assembled constraint equations, corresponding to the rows of C.
    Index dofs{};           ///< Number of active system DOFs, corresponding to the columns of C.
    Index rank{};           ///< Numerically determined rank of the constraint matrix C.
    Index redundant_rows{}; ///< Number of linearly dependent constraint equations.

    // Numerical diagnostics
    Precision rhs_norm{Precision(0)};      ///< Euclidean norm of the constraint right-hand side d.
    Precision residual_norm{Precision(0)}; ///< Norm of the remaining constraint residual ||C u - d||.

    // System properties
    bool homogeneous{true}; ///< True if the constraint right-hand side d is zero.
    bool feasible{true};    ///< True if the constraint equations are mutually consistent.
    bool rank_known{true};  ///< True if rank contains an explicitly determined numerical rank.
};

} // namespace constraint
} // namespace fem
