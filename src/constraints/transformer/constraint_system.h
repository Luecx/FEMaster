/**
 * @file constraint_system.h
 * @brief Declares the assembled sparse constraint system `C u = d`.
 *
 * The types in this file represent the algebraic constraint system generated
 * from symbolic nodal equations. Each retained equation contributes one row to
 * the sparse matrix `C`, while the corresponding prescribed scalar value is
 * stored in the right-hand-side vector `d`.
 *
 * The assembly process maps node-local degrees of freedom to active global
 * system indices. Terms associated with inactive DOFs are omitted. Equations
 * that contain no remaining active terms can either be retained as explicit
 * zero rows or discarded according to `ConstraintOptions`.
 *
 * The assembled system also preserves the source category of every retained
 * row. This metadata allows later stages to distinguish supports, connectors,
 * couplings, ties and other equation origins when computing reactions or
 * producing diagnostics.
 *
 * @see constraint_system.cpp
 * @see constraint_transformer.h
 * @see ../types/equation.h
 *
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
 * @brief Configures the conversion of symbolic equations into the global constraint system.
 *
 * The options control how equations are treated after inactive nodal DOFs have
 * been removed. They affect only the sparse assembly process and do not modify
 * the original symbolic equation collection.
 */
struct ConstraintOptions {
    // Controls whether equations that contain no active DOF entries are retained
    // as explicit zero rows. A positive value enables removal of such rows,
    // while a non-positive value keeps them in the assembled system.
    //
    // The current implementation uses only the sign of this value. Its
    // magnitude is not applied as a numerical coefficient tolerance.
    Precision zero_row_drop_tolerance{Precision(0)};
};

/**
 * @brief Stores the assembled global linear constraint system.
 *
 * The system represents the equations
 *
 *     C u = d,
 *
 * where `u` contains all active global system DOFs. The sparse matrix `C` has
 * one row per retained constraint equation and one column per active DOF.
 *
 * The ordering of `d` and `row_sources` matches the row ordering of `C`
 * exactly. This invariant allows every algebraic row to be traced back to its
 * symbolic equation category.
 */
struct ConstraintSystem {
    // Number of retained constraint equations. This is the row count of C, the
    // length of d and the required number of entries in row_sources.
    Index equations{};

    // Number of active global system DOFs addressed by the constraint matrix.
    // This is the column count of C.
    Index dofs{};

    // Sparse global constraint matrix. Entry C(i, j) contains the coefficient
    // of active global DOF j in retained constraint equation i.
    SparseMatrix C{};

    // Right-hand-side vector of prescribed constraint values. Entry d(i)
    // belongs to row i of C.
    DynamicVector d{};

    // Source category of each retained constraint row. The ordering matches C
    // and d and is used for diagnostics and source-specific force recovery.
    std::vector<EquationSourceKind> row_sources{};
};

// Assemble symbolic nodal equations into the sparse global system
//
//     C u = d.
//
// Every equation entry is mapped through `system_dofs` from a node-local DOF to
// an active global column index. Entries mapped to negative indices are treated
// as inactive and omitted. Retained equations are packed consecutively, so
// their assembled row indices can differ from their original positions when
// empty rows are discarded.
//
// The returned object contains the sparse coefficient matrix, the prescribed
// right-hand side and the source metadata associated with every retained row.
ConstraintSystem assemble_constraint_system(
    const Equations&         equations,
    const SystemDofIds&      system_dofs,
    Index                    n_dofs,
    const ConstraintOptions& options = {}
);

} // namespace constraint
} // namespace fem

