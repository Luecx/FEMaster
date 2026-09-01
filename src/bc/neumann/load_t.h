/**
 * @file load_t.h
 * @brief Declares structural thermal loading driven by a nodal temperature field.
 *
 * `TLoad` does not prescribe equivalent nodal forces directly. Instead, it
 * validates a scalar nodal temperature field and delegates thermal-strain
 * evaluation and force assembly to every structural element, together with the
 * reference temperature that defines the stress-free state.
 *
 * Despite consuming a thermal field, `TLoad` belongs to the structural Neumann
 * hierarchy because its result is an equivalent mechanical right-hand-side
 * contribution. It is distinct from prescribed-temperature Dirichlet conditions
 * used by a thermal conduction analysis.
 *
 * @see TLoad
 * @see StructuralNeumann
 * @see model::StructuralElement
 * @see load_t.cpp
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

/**
 * @brief Applies element-specific thermal expansion loads to the structural model.
 *
 * The temperature input must be a one-component nodal field. Each structural
 * element gathers the temperatures of its own nodes, computes the temperature
 * difference relative to `ref_temp_` and assembles the corresponding equivalent
 * nodal force according to its material and kinematic formulation.
 */
struct TLoad : StructuralNeumann {
    // Shared ownership type for direct references to structural thermal loads.
    using Ptr = std::shared_ptr<TLoad>;

    // Scalar nodal temperature field consumed by structural elements during
    // thermal-expansion load assembly.
    SPtr<model::Field> temp_field_ = nullptr;

    // Temperature at which the structural element is considered free of thermal
    // strain and therefore carries no thermal-expansion stress.
    Precision ref_temp_ = NAN;

    // Validate `temp_field_` and invoke element-specific thermal-load assembly on
    // every structural element. Time and amplitude scaling are intentionally not
    // used because the complete temperature state is represented by the field.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    // Return the temperature-field name and stress-free reference temperature.
    std::string str() const override;
};

} // namespace fem::bc
