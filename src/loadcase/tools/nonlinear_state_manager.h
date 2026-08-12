/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material and contact-state management.
 *
 * Material history is committed per accepted physical increment. Node-to-surface
 * contact additionally owns discrete slave-to-master connectivity. That topology
 * is updated once per Newton residual/tangent evaluation, frozen for nested line
 * search evaluations and checked after convergence for a discontinuity update.
 *
 * Following CalculiX N2F contact, contact-element regeneration is disabled after
 * eight completed Newton iterations of one increment attempt. Continuous contact
 * geometry on the retained master faces is still evaluated at the current nodal
 * positions.
 *
 * @see model::ModelData::material_state
 * @see tools::LoadControl
 * @see tools::ArcLengthControl
 *
 * @author Finn Eggers
 * @date 12.08.2026
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
 * @brief Owns nonlinear material history and coordinates contact topology.
 *
 * Material evaluations overwrite an active trial field in place and therefore
 * restart from committed history for every independent residual or tangent
 * evaluation. Accepted increments promote the trial constitutive state to the
 * committed field.
 *
 * Contact keeps its own transactional partner map. This manager only coordinates
 * when those maps may update, when nested evaluations must freeze them and when
 * the CalculiX iteration-limit freeze becomes active.
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

    // Contact topology lifecycle used by nonlinear path control
    void begin_contact_increment_trial();
    void begin_contact_update_trial();
    void begin_contact_frozen_trial();
    void commit_contact_trial();
    void rollback_contact_trial();
    void finish_contact_iteration();
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

    // CalculiX N2F regenerates contact elements through iteration 8 and freezes
    // their discrete connectivity for later iterations of the same attempt.
    Index contact_iterations_ = 0;
    bool  contact_frozen_     = false;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
