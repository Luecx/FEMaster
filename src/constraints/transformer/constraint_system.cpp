/**
 * @file constraint_system.cpp
 * @brief Implements assembly of the linear constraint system `C u = d`.
 *
 * @see src/constraints/transformer/constraint_system.h
 * @author Finn Eggers
 * @date 14.07.2026
 */

#include "constraint_system.h"

#include "../../core/logging.h"

namespace fem {
namespace constraint {

ConstraintSystem assemble_constraint_system(
    const Equations&         equations,
    const SystemDofIds&      system_dofs,
    Index                    n_dofs,
    const ConstraintOptions& options
) {
    ConstraintSystem system{};
    system.dofs = n_dofs;

    // Prepare sparse matrix and row metadata
    TripletList matrix_entries{};
    matrix_entries.reserve(64 * equations.size());

    std::vector<Precision> rhs_values{};
    rhs_values.reserve(equations.size());

    system.row_sources.reserve(equations.size());

    // Assemble all equations containing active DOFs
    Index equation_id = 0;

    for (const Equation& equation : equations) {
        bool has_active_entry = false;

        for (const EquationEntry& entry : equation.entries) {
            const int dof_id = system_dofs(entry.node_id, entry.dof);

            if (dof_id < 0) {
                continue;
            }

            matrix_entries.emplace_back(equation_id, dof_id, entry.coeff);
            has_active_entry = true;
        }

        const bool keep_equation =
            has_active_entry ||
            options.zero_row_drop_tolerance <= Precision(0);

        if (!keep_equation) {
            continue;
        }

        system.row_sources.push_back(equation.source);
        rhs_values.push_back(equation.rhs);

        ++equation_id;
    }

    // Build the sparse constraint matrix
    system.equations = equation_id;

    system.C.resize(system.equations, system.dofs);
    system.C.setFromTriplets(matrix_entries.begin(), matrix_entries.end());
    system.C.makeCompressed();

    // Build the constraint right-hand side
    system.d.resize(system.equations);

    for (Index row = 0; row < system.equations; ++row) {
        system.d[row] = rhs_values[static_cast<std::size_t>(row)];
    }

    // Report the assembled system
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
