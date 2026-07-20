/**
 * @file constraint_system.cpp
 * @brief Implements assembly of the global sparse constraint system `C u = d`.
 *
 * The symbolic constraint equations are converted into a sparse matrix acting
 * on the active global system DOFs. Each retained symbolic equation contributes
 * one row to `C`, while its prescribed scalar value is stored in the
 * corresponding entry of `d`.
 *
 * Terms that reference inactive nodal DOFs are ignored during assembly. An
 * equation whose terms are all inactive can either be retained as an explicit
 * zero row or removed, depending on the supplied `ConstraintOptions`.
 *
 * The original equation source is retained for every assembled row so later
 * stages can distinguish supports, connectors, couplings, ties and other
 * constraint categories when recovering reactions or producing diagnostics.
 *
 * @see constraint_system.h
 * @see constraint_transformer.h
 * @see ../types/equation.h
 *
 * @author Finn Eggers
 * @date 14.07.2026
 */

#include "constraint_system.h"

#include "../../core/logging.h"

namespace fem {
namespace constraint {

// Assemble a collection of symbolic nodal equations into the sparse global
// constraint system
//
//     C u = d.
//
// The nodal DOF map converts every equation entry into an active global column
// index. Entries belonging to inactive DOFs are omitted. Retained equations are
// packed consecutively, so the assembled row index may differ from the original
// position in the input collection when empty rows are discarded.
ConstraintSystem assemble_constraint_system(
    const Equations&         equations,
    const SystemDofIds&      system_dofs,
    Index                    n_dofs,
    const ConstraintOptions& options
) {
    // Start with an empty system and record the number of active global DOFs.
    // This becomes the column count of the constraint matrix C.
    ConstraintSystem system{};
    system.dofs = n_dofs;

    // Collect sparse matrix coefficients as triplets before constructing C.
    // The reserve value anticipates equations with several coupled nodal terms
    // and reduces reallocations during assembly.
    TripletList matrix_entries{};
    matrix_entries.reserve(64 * equations.size());

    // Store right-hand-side values temporarily because the final number of
    // retained rows is not known until all equations have been inspected.
    std::vector<Precision> rhs_values{};
    rhs_values.reserve(equations.size());

    // Preserve one source identifier for each retained matrix row. This
    // metadata is later used, for example, to isolate support reactions from
    // other constraint forces.
    system.row_sources.reserve(equations.size());

    // Consecutive row index assigned only to equations that survive the
    // inactive-DOF and zero-row filtering below.
    Index equation_id = 0;

    // Process the symbolic equations in their original order so the retained
    // matrix rows and source metadata remain deterministic.
    for (const Equation& equation : equations) {
        // Track whether at least one term maps to an active global DOF. If no
        // active term remains, the equation would become a zero row in C.
        bool has_active_entry = false;

        // Convert each nodal equation term into a sparse matrix coefficient.
        for (const EquationEntry& entry : equation.entries) {
            // Map the node-local DOF pair to its active global system index.
            // Negative values indicate that this DOF is not part of the
            // current active algebraic system.
            const int dof_id = system_dofs(entry.node_id, entry.dof);

            // Ignore terms associated with inactive DOFs. They do not receive a
            // column in the assembled matrix.
            if (dof_id < 0) {
                continue;
            }

            // Insert the coefficient into the current retained row and the
            // mapped global DOF column.
            matrix_entries.emplace_back(equation_id, dof_id, entry.coeff);
            has_active_entry = true;
        }

        // Retain equations containing at least one active term. Equations that
        // collapse to zero rows are retained only when zero-row dropping is
        // disabled by a non-positive tolerance.
        const bool keep_equation =
            has_active_entry ||
            options.zero_row_drop_tolerance <= Precision(0);

        // Skip the complete equation when it contains no active terms and the
        // selected options request removal of such rows.
        if (!keep_equation) {
            continue;
        }

        // Store the source category in the same order as the assembled matrix
        // rows. This establishes the invariant
        //
        //     row_sources.size() == equations.
        system.row_sources.push_back(equation.source);

        // Store the scalar prescribed value associated with this retained row.
        rhs_values.push_back(equation.rhs);

        // Advance the packed matrix row index only after an equation has been
        // accepted.
        ++equation_id;
    }

    // The number of retained equations defines the row count of C and the
    // length of d.
    system.equations = equation_id;

    // Construct the sparse constraint matrix from the collected triplets.
    // Duplicate triplets, if present, are combined by Eigen during insertion.
    system.C.resize(system.equations, system.dofs);
    system.C.setFromTriplets(matrix_entries.begin(), matrix_entries.end());

    // Convert the matrix into compressed sparse storage for efficient
    // factorization, multiplication and traversal in subsequent algorithms.
    system.C.makeCompressed();

    // Allocate the right-hand-side vector with exactly one entry per retained
    // constraint row.
    system.d.resize(system.equations);

    // Copy the temporarily stored scalar values into Eigen's dynamic vector
    // while preserving the same packed row order used for C.
    for (Index row = 0; row < system.equations; ++row) {
        system.d[row] = rhs_values[static_cast<std::size_t>(row)];
    }

    // Report the dimensions and sparsity of the completed constraint system:
    //
    //     m   = number of equations,
    //     n   = number of active DOFs,
    //     nnz = number of stored nonzero coefficients.
    logging::info(
        true,
        "[ConstraintSystem] Assembled C: m=", system.equations,
        " n=",                              system.dofs,
        " nnz=",                            static_cast<Index>(system.C.nonZeros())
    );

    return system;
}

} // namespace constraint
} // namespace fem
