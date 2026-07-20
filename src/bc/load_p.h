/**
 * @file load_p.h
 * @brief Declares scalar pressure loading on surface regions.
 *
 * `PLoad` applies a uniform scalar pressure normal to every selected surface.
 * The surface geometry supplies the position-dependent normal and performs the
 * consistent integration into nodal forces. An optional amplitude scales the
 * pressure magnitude at the current analysis time.
 *
 * @see PLoad
 * @see Load
 * @see load_p.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "load.h"

#include "../data/region.h"

namespace fem {
namespace bc {

/**
 * @brief Integrates a uniform follower direction based on surface normals.
 *
 * The stored pressure is multiplied by the optional amplitude and applied in
 * the negative direction of the surface normal returned at each integration
 * point. Consequently, the vector direction follows the geometric orientation
 * of each individual surface while the scalar magnitude remains uniform.
 */
struct PLoad : public Load {
    // Shared ownership type for direct references to pressure loads.
    using Ptr = std::shared_ptr<PLoad>;

    // Nominal scalar pressure magnitude. The implementation multiplies this
    // value by the optional amplitude before surface integration.
    Precision pressure_ = NAN;

    // Surface region receiving the pressure. Each valid surface is integrated
    // independently into the common nodal load field.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    // Construct an unconfigured pressure load for later population.
    PLoad() = default;

    // Enable polymorphic destruction through `Load`.
    ~PLoad() override = default;

    // Evaluate the optional amplitude, obtain the local surface normal at every
    // integration point and accumulate the resulting consistent nodal forces.
    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;

    // Describe the target surface set, nominal pressure and optional amplitude.
    std::string str() const override;
};
} // namespace bc
} // namespace fem
