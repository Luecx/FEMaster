/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements nonlinear material and contact-state management.
 *
 * Material history is represented by committed and trial material-point fields.
 * Independent residual and tangent evaluations restart from committed material
 * history and may overwrite the active trial field without modifying the last
 * accepted constitutive state.
 *
 * Node-to-surface contact owns transactional slave-to-master connectivity. The
 * existing nonlinear trial callbacks are forwarded directly to every Contact
 * object so increment cutbacks, line-search trials and post-Newton discontinuity
 * checks use the same transaction structure as the generic nonlinear solver.
 *
 * @see NonlinearStateManager
 * @see LoadControl
 * @see ArcLengthControl
 *
 * @author Finn Eggers
 * @date 12.08.2026
 */

#include "nonlinear_state_manager.h"

#include "../../model/model.h"

#include <utility>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * Creates the nonlinear state coordinator for one load-case execution.
 */
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

/**
 * Restores the material-state field that was visible before this manager was
 * created.
 */
NonlinearStateManager::~NonlinearStateManager() {
    model_._data->material_state = previous_material_state_;
}

/**
 * Restarts the active constitutive trial from the last committed material state.
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
 * Opens a contact transaction in which Newton tangent assembly may regenerate
 * the discrete slave-to-master connectivity. Contact itself suppresses updates
 * when its CalculiX regeneration limit has already been reached.
 */
void NonlinearStateManager::begin_contact_update_trial() {
    for (auto& contact : model_._data->contacts) {
        contact.begin_update_trial();
    }
}

/**
 * Opens a nested frozen transaction. Line-search residuals therefore reuse the
 * exact master-face connectivity selected by the tangent that produced the
 * current Newton direction.
 */
void NonlinearStateManager::begin_contact_frozen_trial() {
    for (auto& contact : model_._data->contacts) {
        contact.begin_frozen_trial();
    }
}

/**
 * Commits the current contact transaction.
 */
void NonlinearStateManager::commit_contact_trial() {
    for (auto& contact : model_._data->contacts) {
        contact.commit_trial();
    }
}

/**
 * Restores the discrete contact state saved at trial entry.
 */
void NonlinearStateManager::rollback_contact_trial() {
    for (auto& contact : model_._data->contacts) {
        contact.rollback_trial();
    }
}

/**
 * Reports the CalculiX-style contact-element count flag produced by the most
 * recent Newton regeneration. A significant count change requests another
 * Newton iteration at the same load factor.
 */
bool NonlinearStateManager::update_contact_active_set() {
    bool unchanged = true;
    for (const auto& contact : model_._data->contacts) {
        unchanged = unchanged && !contact.contact_topology_changed();
    }
    return unchanged;
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
