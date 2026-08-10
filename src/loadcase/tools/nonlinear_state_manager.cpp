/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements nonlinear material and contact state management.
 *
 * Material history is represented by committed and trial material-point fields.
 * Contact history remains inside each contact object and is coordinated here so
 * nonlinear load cases do not need to traverse or manipulate contact state
 * directly.
 *
 * @see NonlinearStateManager
 * @see constraint::Contact
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
        contact.reset_augmented_lagrange_penalty();
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

/**
 * Opens the outer contact trial for one attempted nonlinear increment or one
 * post-convergence contact refresh.
 *
 * NodeRegion contact may select one new discrete master partner before freezing
 * it. SurfaceRegion mortar only copies transactional state; its overlap geometry
 * is recomputed on every assembly and is never frozen.
 */
void NonlinearStateManager::begin_contact_update_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.begin_update_trial();
    }
}

/**
 * Opens a nested predictor/line-search trial.
 *
 * NodeRegion partner ownership is frozen in the nested trial. SurfaceRegion
 * mortar again uses the trial only for commit/rollback and recomputes overlap at
 * the candidate geometry.
 */
void NonlinearStateManager::begin_contact_frozen_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.begin_frozen_trial();
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
 * Updates contact state after a converged inner Newton solve.
 *
 * NodeRegion contact still owns a discrete master-facet state. A changed partner
 * signature is therefore re-equilibrated before its multiplier update. Dual
 * mortar contact has no discrete partner state: its nodal unilateral active set
 * is part of the semismooth Newton law, and the outer update acts only on the
 * converged nodal mortar multipliers.
 *
 * Penalty adaptation is common to both formulations and is based on the gaps
 * stored by the corresponding discrete contact constraints.
 *
 * @return `true` when no discrete node-contact partner and no AL multiplier
 *         changed, otherwise `false` so the path controller repeats Newton.
 */
bool NonlinearStateManager::update_contact_active_set() {
    bool changed = false;

    for (const auto& contact : model_._data->contacts) {
        // Mortar coupling has no single master partner to freeze or refresh.
        const bool partner_changed =
            !contact.uses_slave_surface() && contact.partner_signature_changed();

        if (!partner_changed && contact.adapt_augmented_lagrange_penalty()) {
            logging::info(
                true,
                "CONTACT: penetration stagnated; increasing effective penalty to ",
                contact.current_penalty()
            );
        }

        const bool multiplier_changed =
            contact.uses_slave_surface()
                ? contact.update_augmented_lagrange_surface()
                : contact.update_augmented_lagrange();

        changed = changed || partner_changed || multiplier_changed;
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
