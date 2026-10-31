/**
 * @file constraint_groups.cpp
 * @brief Implements utilities for grouped constraint equations.
 *
 * Constraint groups preserve the origin of equations collected from supports,
 * connectors, MPCs and other model constraints. This file provides the ordered
 * equation list needed for assembly, a nodal support mask for result output,
 * and a diagnostic overview of the original groups.
 *
 * @see src/constraints/constraint_groups.h
 */

#include "constraint_groups.h"

#include "transformer/constraint_report.h"
#include "../core/logging.h"

namespace fem {
namespace constraint {

/**
 * @brief Combines all constraint groups into a single ordered equation list.
 *
 * The solver consumes a flat sequence of equations rather than separate
 * categories. The insertion order follows the group order defined by the
 * model assembly and is preserved so that row ordering and source metadata
 * remain consistent with the original constraints.
 *
 * @return A new vector containing all equations in group order.
 */
Equations ConstraintGroups::flatten() const {
    Equations all;

    // Append each group in the same order used during constraint assembly.
    // Keep the groups themselves unchanged for diagnostics and result recovery.
    all.insert(all.end(), supports  .begin(), supports  .end());
    all.insert(all.end(), connectors.begin(), connectors.end());
    all.insert(all.end(), mpcs      .begin(), mpcs      .end());
    all.insert(all.end(), couplings .begin(), couplings .end());
    all.insert(all.end(), ties      .begin(), ties      .end());
    all.insert(all.end(), rbms      .begin(), rbms      .end());
    all.insert(all.end(), pretensions.begin(), pretensions.end());
    all.insert(all.end(), others    .begin(), others    .end());

    return all;
}

/**
 * @brief Creates a nodal mask for degrees of freedom referenced by supports.
 *
 * Support equations may involve several nodal degrees of freedom. Every active
 * entry appearing in such an equation is marked in the mask, while inactive
 * DOFs and entries belonging exclusively to other constraint groups remain
 * unmarked. The mask has the same dimensions as the system DOF index matrix
 * and is used to select support reaction components for result output.
 *
 * @param dofs Global nodal DOF indices; -1 denotes an inactive DOF.
 * @return Boolean matrix marking active DOFs referenced by support equations.
 */
BooleanMatrix ConstraintGroups::support_mask(const SystemDofIds& dofs) const {
    BooleanMatrix mask = BooleanMatrix::Constant(dofs.rows(), dofs.cols(), false);

    // Traverse only support equations, excluding other constraint families.
    for (const auto& equation : supports) {
        for (const auto& entry : equation.entries) {
            // Ignore out-of-range nodes/DOFs and inactive system coordinates.
            if (entry.node_id >= 0 && entry.node_id < dofs.rows()
                && entry.dof < dofs.cols()
                && dofs(entry.node_id, entry.dof) != -1) {
                mask(entry.node_id, entry.dof) = true;
            }
        }
    }

    return mask;
}

/**
 * @brief Logs the collected equations grouped by their physical origin.
 *
 * The overview shows the number of equations for every category, including
 * empty groups. For nonempty groups, each equation is formatted using the
 * common constraint formatting utility and printed with nested indentation.
 * The stored equations and their ordering are not modified.
 */
void ConstraintGroups::report() const {
    using fem::logging::info;

    const EquationFormatOptions format{};

    // Associate each group with its display label in the standard report order.
    const struct {
        const char* name;
        const Equations& equations;
    } categories[] = {
        {"Supports",   supports},
        {"Connectors", connectors},
        {"MPCs",       mpcs},
        {"Couplings",  couplings},
        {"Ties",       ties},
        {"RBMs",       rbms},
        {"Pretensions", pretensions},
        {"Other",      others}
    };

    // Always show category counts; print individual equations only when present.
    info(true, "Constraint overview:");
    logging::up();
    for (const auto& category : categories) {
        info(true, category.name, " (", category.equations.size(), ")");
        if (!category.equations.empty()) {
            logging::up();
            for (const auto& line : format_equations(category.equations, format)) {
                info(true, "    ", line);
            }
            logging::down();
        }
    }
    logging::down();
}

} // namespace constraint
} // namespace fem
