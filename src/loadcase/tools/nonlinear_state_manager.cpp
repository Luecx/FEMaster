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

NonlinearStateManager::NonlinearStateManager(model::Model& model)
    : model_(model),
      previous_material_state_(model._data->material_state) {
    for (const auto& contact : model_._data->contacts) {
        contact.reset_penalty_adaptation();
    }

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

void NonlinearStateManager::begin_contact_trial() {
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

void NonlinearStateManager::bind_material_state() {
    committed_material_state_->name = "MATERIAL_STATE_COMMITTED";
    trial_material_state_->name     = "MATERIAL_STATE";
    model_._data->material_state    = trial_material_state_;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
