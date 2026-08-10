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
 * Opens the outer contact trial of one attempted nonlinear increment or one
 * post-convergence multiplier refresh.
 *
 * Surface mortar has no discrete partner state. The operation therefore only
 * copies the current augmented-Lagrange state transactionally.
 */
void NonlinearStateManager::begin_contact_update_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.begin_trial();
    }
}

/**
 * Opens one nested contact trial for predictor or line-search evaluation.
 *
 * This is intentionally identical to the outer contact trial operation because
 * surface mortar always recomputes geometry from the active nodal positions.
 */
void NonlinearStateManager::begin_contact_frozen_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.begin_trial();
    }
}

void NonlinearStateManager::commit_contact_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.commit_trial();
    }
}

void NonlinearStateManager::rollback_contact_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.rollback_trial();
    }
}

/**
 * Updates augmented-Lagrange multipliers after a converged inner Newton solve.
 *
 * @return `true` when no multiplier changed, otherwise `false` so the path
 *         controller repeats Newton at the same load factor.
 */
bool NonlinearStateManager::update_contact_active_set() {
    bool changed = false;

    for (const auto& contact : model_._data->contacts) {
        changed = contact.update_augmented_lagrange() || changed;
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
