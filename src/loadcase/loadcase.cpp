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
 * and contact multipliers. An independent analysis restores reference geometry
 * and discards material history, which is lazily recreated before element
 * evaluation. Explicit initial temperature and velocity remain unchanged.
 */
void LoadCase::prepare_state() const {
    if (resume) return;

    *model->_data->positions = *model->_data->positions_reference;
    const auto& reference_temperature = model->_data->temperature_reference;
    model->_data->temperature = reference_temperature
        ? std::make_shared<model::Field>(*reference_temperature) : nullptr;
    model->_data->material_state_old.reset();
    model->_data->material_state_new.reset();

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
 * Assembles the static transition between initial and target loads.
 *
 * Shared definitions retain their nominal magnitude. Removed definitions ramp
 * down linearly; added definitions follow their step-time amplitude, or ramp
 * linearly by default. The same assembly handles all structural load families.
 */
model::Field LoadCase::build_step_load_matrix(Precision time, Precision progress) const {
    return model->build_load_matrix(
        initial_conditions, target_conditions,
        time, progress);
}

/**
 * Evaluates only the target definitions with time-dependent amplitudes.
 *
 * Transient and frequency-domain analyses do not perform static load ramps.
 * With an empty initial state, the target conditions use their time-dependent
 * amplitude or the unit progress factor.
 */
model::Field LoadCase::build_target_load_matrix(Precision time) const {
    return model->build_load_matrix(
        bc::ConditionState{}, target_conditions,
        time, Precision(1));
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
