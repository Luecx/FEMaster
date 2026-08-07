/**
 * @file nonlinear_state_manager.cpp
 * @brief Implements transactional nonlinear state management.
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
    prepare_evaluation();
}

NonlinearStateManager::~NonlinearStateManager() {
    model_._data->material_state = previous_material_state_;
}

void NonlinearStateManager::prepare_evaluation() {
    if (!committed_material_state_) {
        return;
    }

    trial_material_state_->values = committed_material_state_->values;
    bind_material_state();
}

void NonlinearStateManager::begin_increment_trial() {
    prepare_evaluation();

    for (const auto& contact : model_._data->contacts) {
        contact.begin_increment_trial();
    }
}

void NonlinearStateManager::commit_increment_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.commit_increment_trial();
    }

    if (!committed_material_state_) {
        return;
    }

    std::swap(committed_material_state_, trial_material_state_);
    trial_material_state_->values = committed_material_state_->values;
    bind_material_state();
}

void NonlinearStateManager::rollback_increment_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.rollback_increment_trial();
    }

    prepare_evaluation();
}

void NonlinearStateManager::begin_line_search_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.begin_line_search_trial();
    }
}

void NonlinearStateManager::commit_line_search_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.commit_line_search_trial();
    }
}

void NonlinearStateManager::rollback_line_search_trial() {
    for (const auto& contact : model_._data->contacts) {
        contact.rollback_line_search_trial();
    }
}

bool NonlinearStateManager::update_active_set(const Evaluation& evaluation) {
    if (model_._data->contacts.empty()) {
        return true;
    }

    for (const auto& contact : model_._data->contacts) {
        contact.begin_active_set_trial();
    }

    try {
        evaluation();
    } catch (...) {
        for (const auto& contact : model_._data->contacts) {
            contact.rollback_active_set_trial();
        }
        throw;
    }

    bool changed = false;

    for (const auto& contact : model_._data->contacts) {
        const bool partner_changed    = contact.partner_signature_changed();
        const bool multiplier_changed = contact.update_augmented_lagrange();

        changed = changed || partner_changed || multiplier_changed;
        contact.commit_active_set_trial();
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
