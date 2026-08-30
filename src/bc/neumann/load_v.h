/**
 * @file load_v.h
 * @brief Declares distributed body-force loading over element regions.
 *
 * `VLoad` defines a vector body force per unit mass. Structural elements
 * integrate the field with material density and distribute the resulting force
 * consistently to their nodes. Optional orientation and amplitude modifiers are
 * inherited from `Neumann`.
 *
 * @see VLoad
 * @see Neumann
 * @see load_v.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"

#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates density-scaled body forces over selected elements.
 */
struct VLoad : Neumann {
    using Ptr = std::shared_ptr<VLoad>;

    // Nominal body-force components per unit mass.
    Vec3 values_ = {NAN, NAN, NAN};

    // Element region receiving the distributed load.
    SPtr<model::ElementRegion> region_ = nullptr;

    VLoad() = default;
    ~VLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;
    std::string str() const override;
};

} // namespace fem::bc
