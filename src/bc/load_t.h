/**
 * @file load_t.h
 * @brief Declares thermal loading driven by a nodal temperature field.
 *
 * `TLoad` does not prescribe equivalent nodal forces directly. Instead, it
 * validates a scalar nodal temperature field and delegates thermal-strain
 * evaluation and force assembly to every structural element, together with the
 * reference temperature that defines the stress-free state.
 *
 * @see TLoad
 * @see Load
 * @see load_t.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "load.h"

namespace fem {
namespace bc {

/**
 * @brief Applies element-specific thermal expansion loads to the model.
 *
 * The temperature input must be a one-component nodal field. Each structural
 * element gathers the temperatures of its own nodes, computes the temperature
 * difference relative to `ref_temp_` and assembles the corresponding equivalent
 * nodal force according to its material and kinematic formulation.
 */
struct TLoad : public Load {
    // Shared ownership type for direct references to thermal loads.
    using Ptr = std::shared_ptr<TLoad>;

    // Scalar nodal temperature field consumed by structural elements during
    // thermal-load assembly.
    SPtr<model::Field> temp_field_ = nullptr;

    // Temperature at which the element is considered free of thermal strain.
    Precision ref_temp_ = NAN;

    // Construct an unconfigured thermal load for later collector assignment.
    TLoad() = default;

    // Enable polymorphic destruction through `Load`.
    ~TLoad() override = default;

    // Validate `temp_field_` and invoke `apply_tload()` on every structural
    // element. The current implementation does not use time or amplitude
    // scaling because the complete temperature history is represented by the
    // supplied field.
    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;

    // Return the temperature-field name and reference temperature used by this
    // load definition.
    std::string str() const override;
};
} // namespace bc
} // namespace fem
