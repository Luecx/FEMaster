/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material and contact-state management.
 *
 * Material history is committed per accepted physical increment. Node-to-surface
 * contact additionally owns discrete slave-to-master connectivity. That topology
 * is updated by Newton tangent evaluations, frozen for nested line-search
 * evaluations and checked after convergence for a discontinuity update.
 *
 * The Contact object itself counts CalculiX N2F regeneration iterations and
 * freezes its discrete connectivity after the eighth Newton update assembly.
 * This manager only forwards the existing generic nonlinear trial lifecycle.
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
 * @brief Owns nonlinear material history and coordinates contact transactions.
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
    void begin_contact_update_trial();
    void begin_contact_frozen_trial();
    void commit_contact_trial();
    void rollback_contact_trial();
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
