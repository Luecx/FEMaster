/******************************************************************************
 * @file constraint_set.h
 * @brief Assemble the sparse constraint system `C u = d` from high-level equations.
 *
 * -----------------------------------------------------------------------------
 * ## What is a ConstraintSet?
 *
 * A `ConstraintSet` is the concrete, assembled representation of your linear
 * constraints. You hand it a list of symbolic equations (each equation is a
 * sum of DOF-entries with a right-hand side), plus your global DOF mapping.
 * It produces the sparse matrix `C` and vector `d` such that:
 *
 *     C u = d
 *
 * where `u` stacks **all** system DOFs.
 *
 * -----------------------------------------------------------------------------
 * ## ELI5 example
 *
 * Suppose you want:
 *   - `u( node=5, dof=x ) = 0`
 *   - `u(3,x) - 2 * u(7,x) = 0`
 *   - `u(1,x) + u(8,x) = 1.2`
 *
 * After assembly, you’ll get a 3×n sparse `C` with exactly those coefficients in
 * the right columns (DOFs), and `d = [0, 0, 1.2]^T`.
 *
 * -----------------------------------------------------------------------------
 * ## Scaling options (robustness)
 *
 * Numerical factorizations (QR) are happier when columns of `C` are scaled
 * to comparable magnitudes. `ConstraintSet::Options` exposes:
 *
 *  - `scale_columns`  (default true): scale each column j of `C` to unit 2-norm
 *    when that column is nonzero; store the factor in `col_scale[j]`.
 *  - `scale_rows`     (default false): optional left scaling to tame row
 *    magnitudes; the factors end up in `row_scale[i]`.
 *  - `zero_row_drop_tol`: when > 0 we can drop completely empty rows produced
 *    by equations that reference only inactive DOFs (or after scaling).
 *
 * The builder currently consumes the **scaled** `C` and `d`. If you need to
 * report back in original scaling, you can use `col_scale` / `row_scale`.
 *
 * -----------------------------------------------------------------------------
 * ## Bookkeeping
 *
 *  - `kept_row_ids[i]` tells you which original equation (index into your
 *    `equations` vector) produced the i-th assembled row in `C`.
 *  - `m` and `n` are the assembled sizes (rows, columns / #DOFs).
 *  - `dof_map` is a pointer to your system DOF layout and is not owned.
 *
 * -----------------------------------------------------------------------------
 * @see equation.h for the high-level equation representation
 * @see constraint_builder.h/.cpp for turning (C,d) into a null-space map
 * @date    14.09.2025
 * @author  Finn
 ******************************************************************************/

#pragma once

#include "equation.h"
#include "../core/types_eig.h"
#include "../core/types_cls.h"
#include <vector>

namespace fem::constraint {

struct ConstraintSet {
    // ------------------------- Input ----------------------------------------

    /// High-level equations (each row: sum of entries equals rhs).
    Equations           equations;

    /// Maps (node_id, dof) -> global system DOF index (>= 0) or -1 if inactive.
    /// Not owned; must remain valid during assemble().
    const SystemDofIds* dof_map = nullptr;

    /// Assembly/scaling options.
    struct Options {
        bool      scale_columns     = true;  ///< Column 2-norm scaling (improves QR robustness).
        bool      scale_rows        = false; ///< Optional row scaling (∞-norm), applied via diagonal left-mult.
        Precision zero_row_drop_tol = 0;     ///< >0 to drop fully zero rows; 0 keeps them (as explicit zeros).
    } opt;

    // ------------------------- Output (after assemble) -----------------------

    Index          m = 0;          ///< Effective number of assembled constraint rows.
    Index          n = 0;          ///< Number of global DOFs (columns in C).
    SparseMatrix   C;              ///< Assembled constraint matrix (m × n).
    DynamicVector  d;              ///< Assembled right-hand side (size m).

    /// Row mapping: `kept_row_ids[i]` = index of the original `equations` entry
    /// that produced row i in C (after optional dropping of empty rows).
    std::vector<Index> kept_row_ids;

    /// Column scaling factors (size n). =1 if column scaling disabled or column zero.
    DynamicVector  col_scale;

    /// Row scaling factors (size m). =1 if row scaling disabled.
    DynamicVector  row_scale;

    // ------------------------- API ------------------------------------------

    /**
     * @brief Assemble sparse C and d from `equations` and a system DOF map.
     *
     * @param dofs   Global (node_id, dof) → system DOF index map.
     * @param n_dofs Total number of system DOFs (i.e., number of columns in C).
     *
     * ### Behavior
     * - For each Equation:
     *     - For each EquationEntry:
     *         - Look up the global DOF index via `dofs(node_id, dof)`.
     *         - If index >= 0, add a triplet (row, col, coeff) into C.
     *     - Set d[row] = rhs.
     * - If an equation yields no active columns (all entries inactive):
     *     - When `zero_row_drop_tol > 0`, the row is **dropped**.
     *     - Otherwise it is **kept** as an explicit zero row (useful for diagnostics).
     * - If `scale_columns` is true, each nonzero column of C is scaled to unit 2-norm;
     *   factors are stored in `col_scale` and applied in-place to C and d (d unchanged).
     * - If `scale_rows` is true, rows are scaled to unit ∞-norm via a diagonal left-multiply;
     *   the same factors are stored in `row_scale` and applied to both C and d.
     */
    void assemble(const SystemDofIds& dofs, Index n_dofs);
};

} // namespace fem::constraint
