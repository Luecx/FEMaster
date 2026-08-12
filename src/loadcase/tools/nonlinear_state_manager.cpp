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
 * manager coordinates update and frozen trials and applies the CalculiX N2F rule
 * that contact-element regeneration stops after eight completed Newton
 * iterations of one increment attempt.
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

#include <iostream>
#include <utility>

namespace fem {
namespace loadcase {
namespace tools {

namespace {
constexpr Index contact_update_iterations = 8;
}

/**
 * Creates the nonlinear state coordinator for one load-case execution.
 *
 * The caller-visible material-state pointer is saved first and restored by the
 * destructor. If assigned materials require history variables, matching
 * committed and trial `ELEMENT_MP` fields are allocated and initialized in the
 * global material-point enumeration owned by `ModelData`.
 *
 * @param model Model whose nonlinear state is coordinated.
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
 * created. Contact objects remain owned by the model and keep only committed
 * topology after every completed increment transaction.
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
 * Starts one complete increment attempt with contact regeneration enabled.
 *
 * The previous accepted partner maps are snapshotted by each Contact object so a
 * cutback can restore them exactly. The CalculiX iteration counter is local to
 * the attempted increment and therefore restarts from zero after every cutback.
 */
void NonlinearStateManager::begin_contact_increment_trial() {
    contact_iterations_ = 0;
    contact_frozen_     = false;

    for (auto& contact : model_._data->contacts) {
        contact.begin_update_trial();
    }
}

/**
 * Opens a topology-update transaction unless the CalculiX iteration freeze is
 * already active for this increment attempt.
 */
void NonlinearStateManager::begin_contact_update_trial() {
    for (auto& contact : model_._data->contacts) {
        if (contact_frozen_) {
            contact.begin_frozen_trial();
        } else {
            contact.begin_update_trial();
        }
    }
}

/**
 * Opens a nested frozen transaction. Line-search residuals use the master faces
 * selected by the tangent evaluation that produced the current Newton direction.
 */
void NonlinearStateManager::begin_contact_frozen_trial() {
    for (auto& contact : model_._data->contacts) {
        contact.begin_frozen_trial();
    }
}

/**
 * Commits the current contact transaction while restoring its parent update mode.
 */
void NonlinearStateManager::commit_contact_trial() {
    for (auto& contact : model_._data->contacts) {
        contact.commit_trial();
    }
}

/**
 * Restores the complete discrete contact state saved at trial entry.
 */
void NonlinearStateManager::rollback_contact_trial() {
    for (auto& contact : model_._data->contacts) {
        contact.rollback_trial();
    }
}

/**
 * Advances the CalculiX N2F contact-iteration counter.
 *
 * CalculiX regenerates node-to-face contact elements while `iit <= 8` and keeps
 * them fixed for later Newton iterations. This method is called once after every
 * completed Newton iteration, including iterations belonging to a same-load
 * discontinuity restart.
 */
void NonlinearStateManager::finish_contact_iteration() {
    if (model_._data->contacts.empty() || contact_frozen_) {
        return;
    }

    ++contact_iterations_;
    if (contact_iterations_ < contact_update_iterations) {
        return;
    }

    contact_frozen_ = true;
    for (auto& contact : model_._data->contacts) {
        contact.freeze_partner_updates();
    }

    std::cout << "[CONTACT] topology frozen after " << contact_iterations_
              << " Newton iterations\n";
}

/**
 * Reports whether the post-Newton contact search reproduced the discrete
 * topology used by the converged Newton residual.
 *
 * A changed slave-to-master signature requests one more Newton solve at the same
 * load factor. Once the CalculiX iteration freeze is active no further topology
 * regeneration is permitted and the active set is therefore considered fixed.
 */
bool NonlinearStateManager::update_contact_active_set() {
    if (contact_frozen_) {
        return true;
    }

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
