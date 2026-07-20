/**
 * @file load_d.h
 * @brief Declares distributed vector tractions on surface regions.
 *
 * `DLoad` represents a traction vector per unit surface area. Each selected
 * surface performs its own shape-function integration and scatters the
 * resulting consistent nodal forces. The vector may be defined globally or in
 * a position-dependent local coordinate system and may be scaled by an
 * amplitude during assembly.
 *
 * @see DLoad
 * @see Load
 * @see load_d.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "load.h"

#include "../data/region.h"

namespace fem {
namespace bc {

/**
 * @brief Integrates a prescribed traction vector over selected surfaces.
 *
 * `values_` stores the nominal traction components. `NaN` marks an omitted
 * component and is converted to zero before integration. Without an assigned
 * orientation, the vector is constant in global coordinates. With an
 * orientation, the same local vector is transformed at every quadrature point,
 * allowing curvilinear or otherwise spatially varying bases.
 */
struct DLoad : public Load {
    // Shared ownership type for code that needs the concrete distributed-load
    // interface instead of the generic `Load::Ptr` type.
    using Ptr = std::shared_ptr<DLoad>;

    // Nominal surface-traction vector. Unspecified components are represented
    // by `NaN` and contribute zero after sanitization.
    Vec3 values_ = {NAN, NAN, NAN};

    // Surface region over which the traction is integrated and distributed to
    // the participating nodes.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    // Construct an unconfigured distributed load for later assignment by a
    // collector or input parser.
    DLoad() = default;

    // Enable polymorphic destruction through `Load`.
    ~DLoad() override = default;

    // Integrate the scaled traction over every valid surface in `region_` and
    // accumulate consistent nodal forces in `bc`. A local orientation, when
    // present, is evaluated separately at each integration point.
    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;

    // Return a one-line description of the target surface set, nominal traction
    // components and optional orientation and amplitude.
    std::string str() const override;
};
} // namespace bc
} // namespace fem
