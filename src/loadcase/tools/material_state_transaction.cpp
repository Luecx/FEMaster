/**
 * @file material_state_transaction.cpp
 * @brief Implements nonlinear material-state buffer management.
 *
 * @see material_state_transaction.h
 */

#include "material_state_transaction.h"

#include "../../model/model.h"

#include <utility>

namespace fem {
namespace loadcase {
namespace tools {

MaterialStateTransaction::MaterialStateTransaction(model::Model& model)
    : model_(model),
      previous_old_(model._data->material_state_old),
      previous_new_(model._data->material_state_new) {
    const Index state_size = model_.maximum_material_state_size();
    if (state_size == 0) {
        return;
    }

    committed_ = model_._data->create_field(
        "MATERIAL_STATE_COMMITTED",
        model::FieldDomain::ELEMENT_MP,
        state_size,
        false,
        false
    );

    trial_ = model_._data->create_field(
        "MATERIAL_STATE_TRIAL",
        model::FieldDomain::ELEMENT_MP,
        state_size,
        false,
        false
    );

    model_.initialize_material_state(*committed_);
    trial_->values = committed_->values;
    bind();
}

MaterialStateTransaction::~MaterialStateTransaction() {
    model_._data->material_state_old = previous_old_;
    model_._data->material_state_new = previous_new_;
}

bool MaterialStateTransaction::active() const {
    return committed_ != nullptr;
}

void MaterialStateTransaction::begin_evaluation() {
    if (!active()) {
        return;
    }

    trial_->values = committed_->values;
    bind();
}

void MaterialStateTransaction::commit_increment() {
    if (!active()) {
        return;
    }

    std::swap(committed_, trial_);
    bind();
}

void MaterialStateTransaction::rollback_increment() {
    begin_evaluation();
}

void MaterialStateTransaction::bind() {
    committed_->name = "MATERIAL_STATE_COMMITTED";
    trial_->name     = "MATERIAL_STATE_TRIAL";

    model_._data->material_state_old = committed_;
    model_._data->material_state_new = trial_;
}

} // namespace tools
} // namespace loadcase
} // namespace fem
