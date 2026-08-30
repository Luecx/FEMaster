/**
 * @file load_d.h
 * @brief Declares distributed vector tractions on surface regions.
 *
 * `DLoad` represents a traction vector per unit surface area. Each selected
 * surface performs shape-function integration and scatters the consistent nodal
 * force. The vector may use a local coordinate system and an amplitude.
 *
 * @see DLoad
 * @see Neumann
 * @see load_d.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"

#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates a prescribed traction vector over selected surfaces.
 */
struct DLoad : Neumann {
    using Ptr = std::shared_ptr<DLoad>;

    // Nominal surface-traction components. NaN denotes an omitted component.
    Vec3 values_ = {NAN, NAN, NAN};

    // Surface region receiving the traction.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    DLoad() = default;
    ~DLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;
    std::string str() const override;
};

} // namespace fem::bc
