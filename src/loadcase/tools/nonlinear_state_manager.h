/**
 * @file nonlinear_state_manager.h
 * @brief Defines nonlinear material and contact state management.
 *
 * Material history uses two explicit buffers. The committed field is published
 * as the immutable constitutive source and the trial field as the separate
 * integration target. Constitutive evaluations therefore never require copying
 * or restoring history between Newton iterations or line-search trials.
 *
 * Contact retains its independent transactional augmented-Lagrange semantics.
 *
 * @see constraint::Contact
 * @see model::ModelData::material_state
 * @see model::ModelData::material_state_target
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
 * @brief Owns accepted/trial material history and coordinates contact state.
 *
 * For every constitutive integration the committed field is read-only and the
 * complete updated state is written to the trial field. Failed Newton states can
 * therefore be discarded simply by overwriting the trial field on the next
 * evaluation. An accepted physical increment is committed by swapping the two
 * field handles; no model-wide history copy is required.
 */
class NonlinearStateManager {
public:
    explicit NonlinearStateManager(model::Model& model);
    ~NonlinearStateManager();

    NonlinearStateManager(const NonlinearStateManager&)            = delete;
    NonlinearStateManager& operator=(const NonlinearStateManager&) = delete;

    // Re-publish the current source/target pair. Kept as a lightweight operation
    // so existing nonlinear controller transactions need no special-case changes.
    void reset_material_state();

    // Promote the complete target state of a converged increment to the immutable
    // source state by swapping the two manager-owned field handles.
    void commit_material_state();

    void begin_contact_trial();
    void commit_contact_trial();
    void rollback_contact_trial();

    bool update_contact_active_set();

private:
    // Publish committed as source and trial as target through ModelData.
    void bind_material_state();

    model::Model& model_;

    // Restore both caller-visible handles when the nonlinear analysis ends.
    model::Field::Ptr previous_material_state_        = nullptr;
    model::Field::Ptr previous_material_state_target_ = nullptr;

    model::Field::Ptr committed_material_state_ = nullptr;
    model::Field::Ptr trial_material_state_     = nullptr;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
