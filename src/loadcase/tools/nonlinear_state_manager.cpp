/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements nonlinear material and contact state management.
 *
 * Material history is represented by committed and trial material-point fields.
 * Surface-to-surface contact owns its augmented-Lagrange state internally and is
 * coordinated here only through generic nested begin/commit/rollback operations.
 *
 * @see NonlinearStateManager
 * @see constraint::Contact
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#include "nonlinear_state_manager.h"

#include "../../model/model.h"

#include <utility>

namespace fem {
namespace loadcase {
namespace tools {

NonlinearStateManager::NonlinearStateManager(model::Model& model)
    : model_(model),
      previous_material_state_(model._data->material_state) {
    const Index state_size = model_.maximum_material_state_size();
    if (state_size == 0) {
        return;
    }

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

NonlinearStateManager::~NonlinearStateManager() {
    model_._data->material_state = previous_material_state_;
}

void NonlinearStateManager::reset_material_state() {
    if (!committed_material_state_) {
        return;
    }

    trial_material_state_->values = committed_material_state_->values;
    bind_material_state();
}

void NonlinearStateManager::commit_material_state() {
    if (!committed_material_state_) {
        return;
    }

    std::swap(committed_material_state_, trial_material_state_);
    trial_material_state_->values = committed_material_state_->values;
    bind_material_state();
}

/**
 * Opens an update-capable contact trial.
 *
 * A trial opened at depth zero is the outer transaction of a new increment
 * attempt. Penalty continuation is reset there so every accepted increment and
 * every cutback restart begins with the user-specified penalty. Nested update
 * trials belong to augmented-Lagrange refreshes at the same load factor and keep
 * the current continuation level.
 */
void NonlinearStateManager::begin_contact_update_trial() {
    if (contact_trial_depth_ == 0) {
        for (const auto& contact : model_._data->contacts) {
            contact.reset_penalty_continuation();
        }
    }

    for (const auto& contact : model_._data->contacts) {
        contact.begin_trial();
    }

    ++contact_trial_depth_;
}

/**
 * Opens one nested contact trial for predictor or line-search evaluation.
 *
 * Temporary evaluations never modify penalty continuation. They inherit the
 * penalty level of their surrounding increment transaction and are subsequently
 * committed or rolled back together with their multiplier state.
 */
void NonlinearStateManager::begin_contact_frozen_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.begin_trial();
    }

    ++contact_trial_depth_;
}

void NonlinearStateManager::commit_contact_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.commit_trial();
    }

    if (contact_trial_depth_ > 0) {
        --contact_trial_depth_;
    }
}

void NonlinearStateManager::rollback_contact_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.rollback_trial();
    }

    if (contact_trial_depth_ > 0) {
        --contact_trial_depth_;
    }
}

/**
 * Updates augmented-Lagrange multipliers after a converged inner Newton solve.
 *
 * Every third successful multiplier augmentation increases the effective contact
 * penalty by one decade. The continuation is local to the current increment and
 * capped inside Contact at 1e4 times the user-specified starting penalty.
 *
 * @return `true` when no multiplier changed, otherwise `false` so the path
 *         controller repeats Newton at the same load factor.
 */
bool NonlinearStateManager::update_contact_active_set() {
    bool changed = false;

    for (const auto& contact : model_._data->contacts) {
        const bool contact_changed = contact.update_augmented_lagrange();

        if (contact_changed) {
            contact.advance_penalty_continuation();
        }

        changed = contact_changed || changed;
    }

    return !changed;
}

void NonlinearStateManager::bind_material_state() {
    committed_material_state_->name = "MATERIAL_STATE_COMMITTED";
    trial_material_state_->name     = "MATERIAL_STATE";
    model_._data->material_state    = trial_material_state_;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
