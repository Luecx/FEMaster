/**
 * @file nonlinear_state_manager.h
 * @brief Declares transactional nonlinear state management.
 *
 * The manager centralizes solver-facing state transitions for nonlinear
 * subsystems. Material history is stored in one committed and one trial
 * `ELEMENT_MP` field while `ModelData::material_state` exposes only the active
 * field used by element and material evaluations. Contact keeps its own nested
 * runtime state but is driven through the same increment and line-search
 * lifecycle.
 */

#pragma once

#include "../../data/field.h"

#include <functional>

namespace fem {
namespace model {
class Model;
}

namespace loadcase {
namespace tools {

class NonlinearStateManager {
public:
    using Evaluation = std::function<void()>;

    explicit NonlinearStateManager(model::Model& model);
    ~NonlinearStateManager();

    NonlinearStateManager(const NonlinearStateManager&)            = delete;
    NonlinearStateManager& operator=(const NonlinearStateManager&) = delete;

    // Rebuild the active material trial field before a residual or tangent evaluation.
    void prepare_evaluation();

    // Complete attempted-increment transaction.
    void begin_increment_trial();
    void commit_increment_trial();
    void rollback_increment_trial();

    // Temporary line-search transaction nested inside the active increment.
    void begin_line_search_trial();
    void commit_line_search_trial();
    void rollback_line_search_trial();

    // Refresh discontinuous nonlinear state after Newton convergence.
    bool update_active_set(const Evaluation& evaluation);

private:
    void bind_material_state();

    model::Model& model_;

    model::Field::Ptr previous_material_state_ = nullptr;
    model::Field::Ptr committed_material_state_ = nullptr;
    model::Field::Ptr trial_material_state_ = nullptr;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
