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
 * As a structural Neumann condition, `DLoad` contributes only to the global
 * structural right-hand side. Surface geometry owns the integration measure and
 * interpolation, while this class owns the physical traction definition and its
 * optional coordinate-system and amplitude modifiers.
 *
 * @see DLoad
 * @see StructuralNeumann
 * @see SurfaceInterface
 * @see load_d.cpp
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates a prescribed traction vector over selected surfaces.
 *
 * `values_` stores the nominal traction components. `NaN` marks an omitted
 * component and is converted to zero before integration. Without an assigned
 * orientation, the vector is constant in global coordinates. With an
 * orientation, the same local vector is transformed at every quadrature point,
 * allowing curvilinear or otherwise spatially varying bases.
 */
struct DLoad : StructuralNeumann {
    // Shared ownership type for code that needs the concrete distributed-load
    // interface instead of the generic Neumann or Load pointer type.
    using Ptr = std::shared_ptr<DLoad>;

    // Nominal surface-traction vector. Unspecified components are represented by
    // NaN and contribute zero after sanitization.
    Vec3 values_ = {NAN, NAN, NAN};

    // Surface region over which the traction is integrated and distributed to
    // the participating nodes.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    // Integrate the scaled traction over every valid surface in `region_` and
    // accumulate consistent nodal forces in `rhs`. A local orientation, when
    // present, is evaluated separately at each integration point.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    // Return a one-line description of the target surface set, nominal traction
    // components and optional orientation and amplitude.
    std::string str() const override;
};

} // namespace fem::bc
