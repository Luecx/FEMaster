/**
 * @file constraint_system.h
 * @brief Declares the assembled linear constraint system `C u = d`.
 *
 * @see src/constraints/transformer/constraint_system.cpp
 * @author Finn Eggers
 * @date 14.07.2026
 */

#pragma once

#include "../../core/types_cls.h"
#include "../../core/types_eig.h"
#include "../types/equation.h"

#include <vector>

namespace fem {
namespace constraint {

/**
 * @brief Controls the assembly of the global constraint system.
 */
struct ConstraintOptions {
    // Row handling
    Precision zero_row_drop_tolerance{Precision(0)}; ///< Positive values discard equations without active DOF entries.
};

/**
 * @brief Stores the assembled linear constraint system `C u = d`.
 */
struct ConstraintSystem {
    // System dimensions
    Index equations{}; ///< Number of assembled constraint equations, corresponding to the rows of C.
    Index dofs{};      ///< Number of active system DOFs, corresponding to the columns of C.

    // Constraint system
    SparseMatrix  C{}; ///< Global constraint matrix.
    DynamicVector d{}; ///< Constraint right-hand side.

    // Equation metadata
    std::vector<EquationSourceKind> row_sources{}; ///< Source of each assembled constraint equation.
};

/**
 * @brief Assembles constraint equations into the sparse system `C u = d`.
 *
 * Terms belonging to inactive DOFs are ignored. Depending on the supplied
 * options, equations without remaining active terms are either retained as
 * zero rows or removed from the assembled system.
 *
 * @param equations Constraint equations to assemble.
 * @param system_dofs Mapping from nodal DOFs to active global DOF indices.
 * @param n_dofs Number of active global DOFs.
 * @param options Constraint assembly options.
 *
 * @return Assembled constraint matrix, right-hand side and row metadata.
 */
ConstraintSystem assemble_constraint_system(
    const Equations&         equations,
    const SystemDofIds&      system_dofs,
    Index                    n_dofs,
    const ConstraintOptions& options = {}
);

} // namespace constraint
} // namespace fem

