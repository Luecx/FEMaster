/**
 * @file load_p.h
 * @brief Declares scalar pressure loading on surface regions.
 *
 * `PLoad` applies a uniform scalar pressure normal to every selected surface.
 * The surface geometry supplies the position-dependent normal and performs the
 * consistent integration into nodal forces. An optional amplitude scales the
 * pressure magnitude at the current analysis time.
 *
 * Pressure is a structural Neumann condition: it contributes only equivalent
 * nodal forces to the structural right-hand side. The sign convention is
 * implemented in `load_p.cpp`, where positive pressure acts opposite to the
 * oriented surface normal returned by the geometric surface implementation.
 *
 * @see PLoad
 * @see StructuralNeumann
 * @see SurfaceInterface
 * @see load_p.cpp
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates a uniform pressure along the geometric surface normal.
 *
 * The stored pressure is multiplied by the optional amplitude and applied in
 * the negative direction of the surface normal returned at each integration
 * point. Consequently, the vector direction follows the geometric orientation
 * of each individual surface while the scalar magnitude remains uniform.
 */
struct PLoad : StructuralNeumann {
    // Shared ownership type for direct references to pressure loads.
    using Ptr = std::shared_ptr<PLoad>;

    // Nominal scalar pressure magnitude. The implementation multiplies this
    // value by the optional amplitude before surface integration.
    Precision pressure_ = NAN;

    // Surface region receiving the pressure. Each valid surface is integrated
    // independently into the common structural nodal right-hand-side field.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    // Evaluate the optional amplitude, obtain the surface normal at every
    // integration point and accumulate the resulting consistent nodal forces.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    // Describe the target surface set, nominal pressure and optional amplitude.
    std::string str() const override;
};

} // namespace fem::bc
