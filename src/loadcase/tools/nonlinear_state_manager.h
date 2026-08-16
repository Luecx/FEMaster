/**
 * @file nonlinear_state_manager.h
 * @brief Owns committed/trial nonlinear material and contact state.
 *
 * Material updates read `committed_material_state_` and write
 * `trial_material_state_`. No committed-to-trial copy is required before an
 * evaluation; every constitutive update must overwrite its complete trial row.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#pragma once

#include "../../model/model.h"

#include <memory>

namespace fem::loadcase::tools {

class NonlinearStateManager {
public:
    explicit NonlinearStateManager(model::Model& model);
    ~NonlinearStateManager();

    NonlinearStateManager(const NonlinearStateManager&) = delete;
    NonlinearStateManager& operator=(const NonlinearStateManager&) = delete;

    void reset_material_state();
    void commit_material_state();

    void begin_contact_update_trial();
    void begin_contact_frozen_trial();
    void commit_contact_trial();
    void rollback_contact_trial();
    bool update_contact_active_set();

private:
    void bind_material_state();

    model::Model& model_;
    model::Field::Ptr previous_material_state_old_;
    model::Field::Ptr previous_material_state_new_;
    model::Field::Ptr previous_material_state_;
    model::Field::Ptr committed_material_state_;
    model::Field::Ptr trial_material_state_;
};

} // namespace fem::loadcase::tools
