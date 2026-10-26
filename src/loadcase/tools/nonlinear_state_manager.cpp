/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements nonlinear material and contact state management.
 *
 * Material history is represented by committed input and trial output fields.
 * Every constitutive evaluation reads only the committed material history and
 * may overwrite the trial field without modifying the last accepted state.
 *
 * Surface-to-surface contact owns its augmented-Lagrange multiplier history and
 * nested trial stack internally. This manager only coordinates begin, commit and
 * rollback operations with the nonlinear path controller. Contact geometry is
 * reconstructed from current nodal positions on every assembly. Adapted
 * contact penalties persist across resumed steps and reset on independent ones.
 *
 * @see NonlinearStateManager
 * @see constraint::Contact
 * @see LoadControl
 * @see ArcLengthControl
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#include "nonlinear_state_manager.h"

#include "../../core/logging.h"
#include "../../model/model.h"

#include <utility>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * Creates the nonlinear state coordinator for one load-case execution.
 *
 * Material history is allocated by Model::step_begin() before element caches
 * are built. This manager coordinates trial and committed state transitions,
 * without resetting history when a new nonlinear analysis starts.
 *
 * @param model Model whose material and contact state is coordinated.
 */
NonlinearStateManager::NonlinearStateManager(model::Model& model)
    : model_(model) {
    material_state_size_ = model_.maximum_material_state_size();
    if (material_state_size_ == 0) {
        return;
    }

    // Reuse the two fields initialized before element step caches are built.
    logging::error(model_._data->material_state_old != nullptr,
        "Nonlinear material input state is not initialized");
    logging::error(model_._data->material_state_new != nullptr,
        "Nonlinear material output state is not initialized");
    logging::error(model_._data->material_state_old != model_._data->material_state_new,
        "Nonlinear material input and output states must use separate storage");

    reset_material_state();
}

/**
 * Restarts the active constitutive trial from the last committed material state.
 *
 * Independent Newton, predictor and line-search evaluations write a fresh trial
 * field while continuing to read the committed state. Copying the committed
 * values also preserves unused components and deterministic trial contents.
 */
void NonlinearStateManager::reset_material_state() {
    if (material_state_size_ == 0) {
        return;
    }

    model_._data->material_state_new->values = model_._data->material_state_old->values;
    name_material_states();
}

/**
 * Accepts the current material trial as the new constitutive history.
 *
 * The committed and trial field objects are swapped instead of copying the
 * accepted data twice. The new trial buffer is then initialized from the accepted
 * state so the next independent evaluation starts from identical history.
 */
void NonlinearStateManager::commit_material_state() {
    if (material_state_size_ == 0) {
        return;
    }

    std::swap(model_._data->material_state_old, model_._data->material_state_new);
    model_._data->material_state_new->values = model_._data->material_state_old->values;
    name_material_states();
}

/**
 * Opens one nested transactional contact trial.
 *
 * Surface mortar has no discrete partner state to freeze or refresh. Predictor,
 * line-search, increment and post-Newton transactions therefore share the same
 * begin operation and differ only in how the nonlinear controller commits or
 * rolls them back.
 */
void NonlinearStateManager::begin_contact_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.begin_trial();
    }
}

/**
 * Commits the innermost contact trial into its parent transaction or committed
 * contact history. Geometry is not involved because it is never stored.
 */
void NonlinearStateManager::commit_contact_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.commit_trial();
    }
}

/**
 * Discards the innermost contact trial and restores the parent multiplier and
 * accepted evaluation state. Current geometry requires no explicit restoration.
 */
void NonlinearStateManager::rollback_contact_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.rollback_trial();
    }
}

/**
 * Updates contact after a converged inner Newton solve.
 *
 * Penalty adaptation is evaluated first from the penetration reached after the
 * previous multiplier augmentation. If penetration did not decrease by at least
 * 20 %, the effective penalty is increased by one decade. The following
 * multiplier update then uses this adapted penalty. The effective penalty was
 * reset only when this manager was constructed and therefore persists across
 * accepted increments.
 *
 * The path controller interprets `false` as a request to repeat Newton at the
 * same load factor with the updated contact history.
 *
 * @return `true` when no contact multiplier changed; `false` when at least one
 *         contact requires another Newton equilibrium solve.
 */
bool NonlinearStateManager::update_contact_active_set() {
    bool changed = false;

    for (const auto& contact : model_._data->contacts) {
        if (contact.adapt_penalty()) {
            logging::info(true,
                "CONTACT: penetration stagnated; increasing effective penalty to ",
                contact.current_penalty());
        }

        const bool contact_changed = contact.update_augmented_lagrange();
        contact.finish_augmentation(contact_changed);
        changed = contact_changed || changed;
    }

    return !changed;
}

/**
 * Restores the semantic names of the committed input and trial output fields
 * after their pointers have been swapped by a material-state commit.
 */
void NonlinearStateManager::name_material_states() {
    model_._data->material_state_old->name = "MATERIAL_STATE_OLD";
    model_._data->material_state_new->name = "MATERIAL_STATE_NEW";
}

} // namespace tools
} // namespace loadcase
} // namespace fem
