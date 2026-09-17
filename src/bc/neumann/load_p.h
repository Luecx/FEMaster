/**
 * @file load_p.h
 * @brief Defines scalar pressure loads acting normal to finite-element surfaces.
 *
 * `PLoad` converts a scalar pressure into a traction aligned with the current
 * geometric surface normal and delegates consistent finite-element integration
 * to the selected surfaces.
 *
 * @see Neumann
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Applies scalar pressure to a surface region.
 *
 * With positive stored pressure `p` and current unit normal `n`, FEMaster uses
 * the traction convention
 *
 *     t = -p n.
 *
 * The equivalent nodal contribution of one surface element is therefore
 *
 *     f_e = -integral_Gamma_e N^T p n dGamma.
 *
 * The optional amplitude scales the scalar pressure before integration. The
 * direction is always determined geometrically from the surface normal, so a
 * separate coordinate-system orientation is not used by this load type.
 */
struct PLoad : Neumann {
    // Types
    using Ptr = std::shared_ptr<PLoad>;

    // Nominal scalar pressure. Positive pressure acts opposite to the current
    // surface normal according to the traction convention t = -p n.
    Precision pressure_ = NAN;

    // Target finite-element surfaces over which the pressure is integrated
    SPtr<model::SurfaceRegion> region_ = nullptr;

    // Construction
    PLoad() = default;
    ~PLoad() override = default;

    // Integrate -N^T p n over every selected surface and accumulate the
    // equivalent nodal force contribution in the supplied RHS field.
    void apply(model::ModelData& model_data, model::Field& rhs,
               Precision time, bool ignore_amplitude = false) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
