/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material and contact state management.
 *
 * The nonlinear state manager owns the material-point history buffers used by
 * nonlinear constitutive evaluations and coordinates the runtime trial states
 * maintained by contact constraints.
 *
 * Material history uses one committed and one active trial `ELEMENT_MP` field.
 * `ModelData::material_state` always points to the active trial field so element
 * and material implementations remain independent of solver-level ownership.
 * Every independent nonlinear evaluation resets this trial field from the last
 * accepted increment, while only an accepted increment promotes it to committed
 * history.
 *
 * Contact keeps its discrete partner, active-set and augmented-Lagrange state in
 * the contact implementation itself. The manager only applies the corresponding
 * update or frozen trial operation consistently to every contact in the model.
 *
 * @see constraint::Contact
 * @see model::ModelData::material_state
 * @see tools::LoadControl
 * @see tools::ArcLengthControl
 *
 * @author Finn Eggers
 * @date 07.08.2026
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
 * Material and contact state intentionally retain separate lifecycle operations
 * because their numerical semantics differ. Material evaluations overwrite an
 * active trial field in place and therefore restart from committed history for
 * every residual or tangent evaluation. Contact instead requires nested trial
 * states so temporary predictor and line-search evaluations can be accepted or
 * discarded independently.
 *
 * The manager owns only the material buffers. Contact state remains owned by the
 * individual `Contact` objects and is accessed through their trial API.
 */
class NonlinearStateManager {
public:
    // Construction binds solver-owned trial storage when at least one assigned
    // material requires history; destruction restores the previous model binding.
    explicit NonlinearStateManager(model::Model& model);
    ~NonlinearStateManager();

    NonlinearStateManager(const NonlinearStateManager&)            = delete;
    NonlinearStateManager& operator=(const NonlinearStateManager&) = delete;

    // Material history used by constitutive evaluations
    void reset_material_state();
    void commit_material_state();

    // Contact trials with either one partner update or immediately frozen partners
    void begin_contact_update_trial();
    void begin_contact_frozen_trial();
    void commit_contact_trial();
    void rollback_contact_trial();

    // Contact active-set and augmented-Lagrange update after Newton convergence
    bool update_contact_active_set();

private:
    // Bind the active material trial field exposed to element and material code
    void bind_material_state();

    // Non-owning model whose contacts and active material-state binding are managed
    model::Model& model_;

    // Material state owned for the lifetime of this nonlinear solution
    model::Field::Ptr previous_material_state_  = nullptr;
    model::Field::Ptr committed_material_state_ = nullptr;
    model::Field::Ptr trial_material_state_     = nullptr;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
