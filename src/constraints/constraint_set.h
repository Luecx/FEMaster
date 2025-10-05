/**
 * @file constraint_set.h
 * @brief Declares the sparse constraint system builder `C u = d`.
 *
 * A `ConstraintSet` collects high-level equations and assembles them into a
 * sparse matrix/vector pair suitable for numerical solvers.
 *
 * @see src/constraints/constraint_set.cpp
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "equation.h"

#include <vector>

namespace fem {
namespace constraint {

/**
 * @struct ConstraintSet
 * @brief Aggregates constraint equations into a sparse representation.
 */
struct ConstraintSet {
    Equations equations;            ///< High-level constraint equations.
    const SystemDofIds* dof_map = nullptr; ///< Pointer to DOF numbering state.

    /**
     * @struct Options
     * @brief Configuration options for constraint assembly.
     */
    struct Options {
        bool scale_columns = false;    ///< Enable column scaling when assembling `C`.
        bool scale_rows = false;       ///< Enable row scaling when assembling `C`.
        Precision zero_row_drop_tol = 0; ///< Drop rows with norm below this tolerance.
    } opt;

    Index m = 0;           ///< Number of assembled constraint equations.
    Index n = 0;           ///< Number of DOFs in the assembled system.
    SparseMatrix C;        ///< Sparse constraint matrix.
    DynamicVector d;       ///< Right-hand side vector.

    std::vector<Index> kept_row_ids; ///< Indices of equations retained after filtering.
    DynamicVector col_scale;         ///< Per-column scaling factors.
    DynamicVector row_scale;         ///< Per-row scaling factors.

    /**
     * @brief Assembles the sparse representation from the stored equations.
     *
     * @param dofs DOF numbering map used to locate matrix columns.
     * @param n_dofs Total number of DOFs in the system.
     */
    void assemble(const SystemDofIds& dofs, Index n_dofs);
};

} // namespace constraint
} // namespace fem
