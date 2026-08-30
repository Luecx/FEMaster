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
 * @see Neumann
 * @see load_p.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"

#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates a uniform pressure as a structural Neumann condition.
 *
 * The stored pressure is multiplied by the optional amplitude and applied in
 * the negative direction of the surface normal returned at each integration
 * point. Consequently, the vector direction follows the geometric orientation
 * of each individual surface while the scalar magnitude remains uniform.
 */
struct PLoad : Neumann {
    using Ptr = std::shared_ptr<PLoad>;

    // Nominal scalar pressure magnitude.
    Precision pressure_ = NAN;

    // Surface region receiving the pressure.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    PLoad() = default;
    ~PLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;
    std::string str() const override;
};

} // namespace fem::bc
