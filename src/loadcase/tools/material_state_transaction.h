/**
 * @file material_state_transaction.h
 * @brief Declares nonlinear material-state buffer management.
 *
 * The transaction owns one committed and one trial `ELEMENT_MP` field for the
 * duration of a nonlinear load case. `ModelData` receives only non-owning shared
 * pointers to the currently bound fields so element and material code can use
 * the same semantic access pattern as other model fields.
 *
 * Every constitutive evaluation reconstructs the trial field from the last
 * committed increment. Accepted line-search candidates remain trial states;
 * only an accepted nonlinear increment swaps the two buffers.
 *
 * @see tools::LoadControl
 * @see tools::ArcLengthControl
 * @see material::Elasticity::evaluate
 */

#pragma once

#include "../../data/field.h"

namespace fem {
namespace model {
struct Model;
}

namespace loadcase {
namespace tools {

class MaterialStateTransaction {
public:
    explicit MaterialStateTransaction(model::Model& model);
    ~MaterialStateTransaction();

    MaterialStateTransaction(const MaterialStateTransaction&)            = delete;
    MaterialStateTransaction& operator=(const MaterialStateTransaction&) = delete;

    [[nodiscard]] bool active() const;

    // Rebuild the writable trial state from the last accepted increment.
    void begin_evaluation();

    // Accept or reject the complete attempted nonlinear increment.
    void commit_increment();
    void rollback_increment();

private:
    void bind();

    model::Model& model_;

    model::Field::Ptr previous_old_ = nullptr;
    model::Field::Ptr previous_new_ = nullptr;

    model::Field::Ptr committed_ = nullptr;
    model::Field::Ptr trial_     = nullptr;
};

} // namespace tools
} // namespace loadcase
} // namespace fem
