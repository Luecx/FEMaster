/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements nonlinear material and contact state management.
 *
 * Material history is represented by committed and trial material-point fields.
 * Independent residual and tangent evaluations restart from committed material
 * history and may overwrite the active trial field without modifying the last
 * accepted constitutive state.
 *
 * Surface-to-surface contact owns its augmented-Lagrange multiplier history and
 * nested trial stack internally. This manager only coordinates begin, commit and
 * rollback operations with the nonlinear path controller. Contact geometry is
 * not part of that history and is reconstructed from current nodal positions on
 * every assembly. Penalty adaptation is initialized once per nonlinear analysis
 * and then persists across accepted load increments.
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
 * The caller-visible material-state pointer is saved first and restored by the
 * destructor. Contact penalty adaptation is reset once here, which returns each
 * contact to its user-provided starting penalty while deliberately allowing later
 * adaptations to persist across load increments and cutbacks within this
 * analysis.
 *
 * If assigned materials require history variables, matching committed and trial
 * `ELEMENT_MP` fields are allocated and initialized in the global material-point
 * enumeration owned by `ModelData`. Models without constitutive history avoid
 * allocating those fields entirely.
 *
 * @param model Model whose material and contact state is coordinated.
 */
NonlinearStateManager::NonlinearStateManager(model::Model& model)
    : model_(model),
      previous_material_state_(model._data->material_state) {
    // Penalty continuation belongs to the complete nonlinear analysis, not to
    // individual increments. Reset it exactly once when the manager is created.
    for (const auto& contact : model_._data->contacts) {
        contact.reset_penalty_adaptation();
    }

    const Index state_size = model_.maximum_material_state_size();
    if (state_size == 0) {
        return;
    }

    // Allocate separate accepted and working constitutive history fields. Both
    // use the same global element/material-point row ordering.
    committed_material_state_ = model_._data->create_field(
        "MATERIAL_STATE_COMMITTED",
        model::FieldDomain::ELEMENT_MP,
        state_size,
        false,
        false
    );

    trial_material_state_ = model_._data->create_field(
        "MATERIAL_STATE",
        model::FieldDomain::ELEMENT_MP,
        state_size,
        false,
        false
    );

    model_.initialize_material_state(*committed_material_state_);
    reset_material_state();
}

/**
 * Restores the material-state field that was visible before this manager was
 * created. The manager-owned committed/trial buffers then leave scope normally.
 */
NonlinearStateManager::~NonlinearStateManager() {
    model_._data->material_state = previous_material_state_;
}

/**
 * Restarts the active constitutive trial from the last committed material state.
 *
 * Independent Newton, predictor and line-search evaluations must not inherit
 * in-place material updates from a previously rejected evaluation. Copying the
 * committed values before rebinding the trial field provides that isolation.
 */
void NonlinearStateManager::reset_material_state() {
    if (!committed_material_state_) {
        return;
    }

    trial_material_state_->values = committed_material_state_->values;
    bind_material_state();
}

/**
 * Accepts the current material trial as the new constitutive history.
 *
 * The committed and trial field objects are swapped instead of copying the
 * accepted data twice. The new trial buffer is then initialized from the accepted
 * state so the next independent evaluation starts from identical history.
 */
void NonlinearStateManager::commit_material_state() {
    if (!committed_material_state_) {
        return;
    }

    std::swap(committed_material_state_, trial_material_state_);
    trial_material_state_->values = committed_material_state_->values;
    bind_material_state();
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
 * Publishes the manager-owned trial material field as the active constitutive
 * history and restores the semantic names of committed and trial buffers.
 */
void NonlinearStateManager::bind_material_state() {
    committed_material_state_->name = "MATERIAL_STATE_COMMITTED";
    trial_material_state_->name     = "MATERIAL_STATE";
    model_._data->material_state    = trial_material_state_;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
