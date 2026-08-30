/**
 * @file neumann.h
 * @brief Declares pure right-hand-side boundary conditions.
 */

#pragma once

#include "../load.h"

#include "../../core/types_num.h"
#include "../../data/field.h"

#include <memory>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Boundary condition contributing only to the global right-hand side.
 */
struct Neumann : Load {
    using Ptr = std::shared_ptr<Neumann>;

    virtual ~Neumann() = default;

    virtual void apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;
};

/** @brief Mechanical Neumann condition assembled into a six-component RHS. */
struct StructuralNeumann : Neumann {
    using Ptr = std::shared_ptr<StructuralNeumann>;
    ~StructuralNeumann() override = default;
};

/** @brief Thermal Neumann condition assembled into a scalar temperature RHS. */
struct ThermalNeumann : Neumann {
    using Ptr = std::shared_ptr<ThermalNeumann>;
    ~ThermalNeumann() override = default;
};

} // namespace fem::bc
