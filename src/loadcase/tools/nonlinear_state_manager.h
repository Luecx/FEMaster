/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material-state management.
 *
 * Material history is committed per accepted physical increment. The current
 * node-to-surface penalty contact formulation is stateless: it reconstructs its
 * complete geometry from current nodal positions on every assembly and owns no
 * multiplier, partner, active-set or transactional history.
 *
 * The contact callback methods remain as no-op adapters for the generic nonlinear
 * path-controller callback interface. They do not store, freeze, switch or update
 * any contact state.
 *
 * @see model::ModelData::material_state
 * @see tools::LoadControl
 * @see tools::ArcLengthControl
 *
 * @author Finn Eggers
 * @date 11.08.2026
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
 * @brief Owns nonlinear material history for one nonlinear analysis.
 *
 * Material evaluations overwrite an active trial field in place and therefore
 * restart from committed history for every independent residual or tangent
 * evaluation. Accepted increments promote the trial constitutive state to the
 * committed field.
 *
 * Node-to-surface penalty contact has no persistent nonlinear state. The contact
 * callback methods are deliberately inert and exist only because the generic
 * nonlinear controller invokes the same lifecycle hooks for any optional
 * stateful subsystem.
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

    // Stateless contact adapters used by generic nonlinear callbacks
    void begin_contact_update_trial();
    void begin_contact_frozen_trial();
    void commit_contact_trial();
    void rollback_contact_trial();
    bool update_contact_active_set();

private:
    // Publish the trial material field as the active ModelData state
    void bind_material_state();

    // Model whose nonlinear material state is coordinated
    model::Model& model_;

    // Caller-visible material state and manager-owned committed/trial buffers
    model::Field::Ptr previous_material_state_  = nullptr;
    model::Field::Ptr committed_material_state_ = nullptr;
    model::Field::Ptr trial_material_state_     = nullptr;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
