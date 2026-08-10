/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material and contact state management.
 *
 * The nonlinear state manager owns material-point history buffers and coordinates
 * the transactional augmented-Lagrange state maintained by surface-to-surface
 * contact. Contact geometry is always reconstructed from current nodal positions
 * and is never frozen or stored between evaluations.
 *
 * Material history is committed per accepted physical increment. Contact
 * multiplier trials follow the nested predictor, line-search and augmentation
 * transactions of the nonlinear controller, while the numerical contact penalty
 * is initialized once for the complete nonlinear analysis and subsequently
 * managed by the contact definition itself.
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
 * evaluation. Accepted increments promote the trial constitutive state to the
 * committed field.
 *
 * Surface contact retains its own nested augmented-Lagrange state. This manager
 * only synchronizes contact begin/commit/rollback operations with the path
 * controller and triggers the post-Newton penalty/multiplier update. No contact
 * geometry or master-facet ownership is stored here.
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

    // Nested contact multiplier/evaluation transactions
    void begin_contact_trial();
    void commit_contact_trial();
    void rollback_contact_trial();

    // Post-Newton augmented-Lagrange multiplier and penalty update
    bool update_contact_active_set();

private:
    // Publish the trial material field as the active ModelData state
    void bind_material_state();

    // Model whose nonlinear state is coordinated
    model::Model& model_;

    // Caller-visible material state and manager-owned committed/trial buffers
    model::Field::Ptr previous_material_state_  = nullptr;
    model::Field::Ptr committed_material_state_ = nullptr;
    model::Field::Ptr trial_material_state_     = nullptr;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
