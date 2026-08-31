/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements nonlinear material and contact state management.
 *
 * Material history is represented by separate committed source and trial target
 * fields. Constitutive integration always reads the committed row and overwrites
 * the corresponding trial row, so rejected residual/tangent evaluations require
 * no material-state rollback or copying.
 *
 * Surface-to-surface contact keeps its own transactional augmented-Lagrange
 * history and remains coordinated independently.
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
      previous_material_state_(model._data->material_state),
      previous_material_state_target_(model._data->material_state_target) {
    // Contact penalty continuation belongs to the complete nonlinear analysis.
    for (const auto& contact : model_._data->contacts) {
        contact.reset_penalty_adaptation();
    }

    const Index state_size = model_.maximum_material_state_size();
    if (state_size == 0) {
        // Stateless models expose no material history buffers at all.
        model_._data->material_state        = nullptr;
        model_._data->material_state_target = nullptr;
        return;
    }

    // Both fields use the same material-point enumeration and width. They are
    // initialized independently once; afterwards integration overwrites the
    // complete target state at every visited material point.
    committed_material_state_ = model_._data->create_field(
        "MATERIAL_STATE_COMMITTED",
        model::FieldDomain::ELEMENT_MP,
        state_size,
        false,
        false
    );

    trial_material_state_ = model_._data->create_field(
        "MATERIAL_STATE_TARGET",
        model::FieldDomain::ELEMENT_MP,
        state_size,
        false,
        false
    );

    model_.initialize_material_state(*committed_material_state_);
    model_.initialize_material_state(*trial_material_state_);
    bind_material_state();
}

NonlinearStateManager::~NonlinearStateManager() {
    model_._data->material_state        = previous_material_state_;
    model_._data->material_state_target = previous_material_state_target_;
}

/**
 * Re-publishes the current source/target pair.
 *
 * No values are copied. Every constitutive integration starts explicitly from
 * the immutable source row and writes a complete target row, so stale values in
 * the uncommitted target field are irrelevant.
 */
void NonlinearStateManager::reset_material_state() {
    if (!committed_material_state_) {
        return;
    }

    bind_material_state();
}

/**
 * Commits the target history of a converged physical increment.
 *
 * Swapping the field handles makes the just-integrated target the new immutable
 * source. The old source becomes scratch target storage for the next increment;
 * it need not be initialized or copied because integration overwrites it.
 */
void NonlinearStateManager::commit_material_state() {
    if (!committed_material_state_) {
        return;
    }

    std::swap(committed_material_state_, trial_material_state_);
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

/**
 * Publishes the explicit source and target material-history fields.
 */
void NonlinearStateManager::bind_material_state() {
    committed_material_state_->name = "MATERIAL_STATE_COMMITTED";
    trial_material_state_->name     = "MATERIAL_STATE_TARGET";

    model_._data->material_state        = committed_material_state_;
    model_._data->material_state_target = trial_material_state_;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
