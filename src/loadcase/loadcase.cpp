/**
 * @file loadcase.cpp
 * @brief Implements the base load-case utilities.
 *
 * Provides diagnostic helpers that report the assembled constraints for a
 * particular load case.
 *
 * @see src/loadcase/loadcase.h
 * @see src/model/constraint_groups.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "loadcase.h"

#include "../constraints/constraint_groups.h"

namespace fem {
namespace loadcase {

/**
 * Prepares the starting configuration of the current analysis.
 *
 * A resumed analysis preserves the prior converged geometry, material history
 * and contact multipliers. Otherwise the reference geometry is restored and
 * material history is marked for initialization by the nonlinear solver.
 * Explicit temperature and initial velocity definitions remain unchanged.
 */
void LoadCase::prepare_state() const {
    if (resume) return;

    *model->_data->positions = *model->_data->positions_reference;
    model->_data->material_state_initialized = false;

    for (const auto& contact : model->_data->contacts) {
        contact.reset_state();
    }
}

/**
 * Stores a converged total nodal displacement in the persistent configuration.
 *
 * Translational coordinates are relative to the immutable reference geometry.
 * Rotational degrees of freedom are stored as total rotation values.
 */
void LoadCase::commit_state(const model::Field& displacement) const {
    auto& current         = *model->_data->positions;
    const auto& reference = *model->_data->positions_reference;

    for (Index node = 0; node < current.rows; ++node) {
        for (Index dof = 0; dof < 3; ++dof) {
            current(node, dof) = reference(node, dof) + displacement(node, dof);
        }
        for (Index dof = 3; dof < 6; ++dof) {
            current(node, dof) = displacement(node, dof);
        }
    }
}

/**
 * @copydoc LoadCase::report_constraint_groups
 */
void LoadCase::report_constraint_groups(const constraint::ConstraintGroups& groups) const {
    if (!report_constraints) {
        return;
    }
    groups.report();
}
} // namespace loadcase
} // namespace fem
