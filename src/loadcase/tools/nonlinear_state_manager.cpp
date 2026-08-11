/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements nonlinear material-state management.
 *
 * Material history is represented by committed and trial material-point fields.
 * Independent residual and tangent evaluations restart from committed material
 * history and may overwrite the active trial field without modifying the last
 * accepted constitutive state.
 *
 * Node-to-surface penalty contact is stateless. The contact lifecycle methods in
 * this file are therefore explicit no-ops used only to satisfy the generic
 * nonlinear path-controller callback interface.
 *
 * @see NonlinearStateManager
 * @see LoadControl
 * @see ArcLengthControl
 *
 * @author Finn Eggers
 * @date 11.08.2026
 */

#include "nonlinear_state_manager.h"

#include "../../model/model.h"

#include <utility>

namespace fem {
namespace loadcase {
namespace tools {

/**
 * Creates the nonlinear material-state coordinator for one load-case execution.
 *
 * The caller-visible material-state pointer is saved first and restored by the
 * destructor. If assigned materials require history variables, matching
 * committed and trial `ELEMENT_MP` fields are allocated and initialized in the
 * global material-point enumeration owned by `ModelData`.
 *
 * @param model Model whose nonlinear material state is coordinated.
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
 * created. The manager-owned committed/trial buffers then leave scope normally.
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
 * Node-to-surface penalty contact owns no update transaction.
 */
void NonlinearStateManager::begin_contact_update_trial() {}

/**
 * Node-to-surface penalty contact owns no frozen transaction.
 */
void NonlinearStateManager::begin_contact_frozen_trial() {}

/**
 * Node-to-surface penalty contact owns no trial state to commit.
 */
void NonlinearStateManager::commit_contact_trial() {}

/**
 * Node-to-surface penalty contact owns no trial state to roll back.
 */
void NonlinearStateManager::rollback_contact_trial() {}

/**
 * Stateless penalty contact has no post-Newton multiplier or active-set update.
 * The current geometry is already reconstructed by every residual/tangent
 * assembly, so the generic controller never needs a same-load Newton restart.
 */
bool NonlinearStateManager::update_contact_active_set() {
    return true;
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
