/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material and contact state management.
 *
 * The nonlinear state manager coordinates the committed and trial material-point
 * history buffers stored by `ModelData` and the transactional augmented-Lagrange
 * state maintained by surface-to-surface contact. Contact geometry is always
 * reconstructed from current nodal positions and is never stored between
 * evaluations.
 *
 * Material history is committed per accepted physical increment. Contact
 * multiplier trials follow the nested predictor, line-search and augmentation
 * transactions of the nonlinear controller, while the numerical contact penalty
 * is initialized once for the complete nonlinear analysis and subsequently
 * managed by the contact definition itself.
 *
 * @see constraint::Contact
 * @see model::ModelData::material_state_old
 * @see model::ModelData::material_state_new
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
struct Model;
}

namespace loadcase {
namespace tools {

/**
 * @brief Coordinates nonlinear material history and contact trial state.
 *
 * Material evaluations read the committed field and write a separate trial
 * field. Independent residual and tangent evaluations therefore cannot advance
 * their own input history. Accepted increments promote the trial constitutive
 * state to the committed field.
 *
 * Surface contact retains its own nested augmented-Lagrange state. This manager
 * only synchronizes contact begin/commit/rollback operations with the path
 * controller and triggers the post-Newton penalty/multiplier update. No contact
 * geometry or master-facet ownership is stored here.
 */
class NonlinearStateManager {
public:
    explicit NonlinearStateManager(model::Model& model);
    ~NonlinearStateManager() = default;

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
    // Restore semantic field names after committed/trial pointer swaps
    void name_material_states();

    // Model whose nonlinear state is coordinated
    model::Model& model_;

    // Uniform constitutive history width; zero denotes stateless materials
    Index material_state_size_ = 0;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
