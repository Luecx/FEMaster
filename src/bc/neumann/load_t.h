/**
 * @file load_t.h
 * @brief Declares structural thermal-expansion loading from a nodal temperature field.
 *
 * `TLoad` does not prescribe equivalent nodal forces directly. It validates a
 * scalar nodal temperature field and delegates thermal-strain evaluation to
 * every structural element together with the stress-free reference temperature.
 *
 * @see TLoad
 * @see Neumann
 * @see load_t.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

/**
 * @brief Applies element-specific thermal expansion loads to the model.
 */
struct TLoad : Neumann {
    using Ptr = std::shared_ptr<TLoad>;

    // Scalar nodal temperature field consumed by structural elements.
    SPtr<model::Field> temp_field_ = nullptr;

    // Stress-free reference temperature.
    Precision ref_temp_ = NAN;

    TLoad() = default;
    ~TLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;
    std::string str() const override;
};

} // namespace fem::bc
