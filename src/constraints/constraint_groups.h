/******************************************************************************
 * @file constraint_groups.h
 * @brief Declares the container that categorises model constraint equations.
 *
 * Groups of constraints are collected during model assembly and can then be
 * reported or flattened into a single list for solver consumption.
 *
 * @see src/constraints/constraint_groups.cpp
 ******************************************************************************/

#pragma once

#include "equation.h"
#include "../core/types_num.h"

#include <vector>

namespace fem {
namespace constraint {

/******************************************************************************
 * @struct ConstraintGroups
 * @brief Stores constraint equations grouped by their origin.
 ******************************************************************************/
struct ConstraintGroups {
    constraint::Equations supports;   ///< Equations originating from supports.
    constraint::Equations connectors; ///< Equations originating from connectors.
    constraint::Equations couplings;  ///< Equations originating from couplings.
    constraint::Equations ties;       ///< Equations originating from tie constraints.
    constraint::Equations others;     ///< Manual or uncategorised constraint equations.

    /// Flattens all groups into a single vector of equations.
    constraint::Equations flatten() const {
        constraint::Equations all;
        auto append = [&all](const constraint::Equations& input) {
            all.insert(all.end(), input.begin(), input.end());
        };
        append(supports);
        append(connectors);
        append(couplings);
        append(ties);
        append(others);
        return all;
    }

    /******************************************************************************
     * @brief Writes a formatted summary of all constraint groups to the logger.
     *
     * @param loadcase_id Identifier of the load case being reported.
     * @param element_dims Spatial dimension of the model's elements.
     ******************************************************************************/
    void report(ID loadcase_id, Dim element_dims) const;
};

} // namespace constraint
} // namespace fem

