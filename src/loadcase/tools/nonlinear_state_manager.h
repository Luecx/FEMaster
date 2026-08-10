/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material and contact state management.
 *
 * The nonlinear state manager owns material-point history buffers and coordinates
 * the transactional augmented-Lagrange state maintained by surface-to-surface
 * contact. Contact geometry is always reconstructed from the current nodal
 * positions and is never frozen or stored between evaluations.
 *
 * @see constraint::Contact
 * @see model::ModelData::material_state
 * @see tools::LoadControl
 * @see tools::ArcLengthControl
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#pragma once

#include "../../data/field.h"

namespace fem {
namespace model {
class Model;
}

namespace loadcase {
namespace tools {

/**
 * @brief Owns nonlinear material history and coordinates contact trial state.
 *
 * Material evaluations overwrite an active trial field in place and therefore
 * restart from committed history for every independent residual or tangent
 * evaluation. Contact owns nested multiplier states so temporary predictor and
 * line-search evaluations can be accepted or discarded independently.
 */
class NonlinearStateManager {
public:
    explicit NonlinearStateManager(model::Model& model);
    ~NonlinearStateManager();

    NonlinearStateManager(const NonlinearStateManager&)            = delete;
    NonlinearStateManager& operator=(const NonlinearStateManager&) = delete;

    // Material history used by constitutive evaluations
    void reset_material_state();
    void commit_material_state();

    // Contact transactions. Update and frozen trials currently have identical
    // state semantics; both names remain because the path controllers distinguish
    // outer active-set refreshes from temporary predictor/line-search evaluations.
    void begin_contact_update_trial();
    void begin_contact_frozen_trial();
    void commit_contact_trial();
    void rollback_contact_trial();

    // Post-Newton augmented-Lagrange multiplier update
    bool update_contact_active_set();

private:
    void bind_material_state();

    model::Model& model_;

    model::Field::Ptr previous_material_state_  = nullptr;
    model::Field::Ptr committed_material_state_ = nullptr;
    model::Field::Ptr trial_material_state_     = nullptr;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
