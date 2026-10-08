/**
 * @file load_p.h
 * @brief Defines scalar pressure loads acting normal to finite-element surfaces.
 *
 * `PLoad` converts a scalar pressure into a traction aligned with the current
 * geometric surface normal and delegates consistent finite-element integration
 * to the selected surfaces.
 *
 * @see Condition
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "../condition.h"
#include "../../constraints/types/equation.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"
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
struct PLoad : Condition {
    // Types
    using Ptr = std::shared_ptr<PLoad>;

    // Target and initial pressure magnitudes per unit surface area.
    // Positive pressure acts opposite to the surface normal (t = -p n).
    // NaN in pressure_start_ indicates an unspecified initial pressure,
    // which is treated as zero during a step transition.
    Precision pressure_       = NAN;
    Precision pressure_start_ = NAN;

    // Target finite-element surfaces over which the pressure is integrated
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Construction
    PLoad() = default;
    ~PLoad() override = default;

    // Integrate the pressure traction over every surface in region_ and
    // accumulate the consistent equivalent nodal forces in the global RHS.
    //
    // Without an amplitude, the pressure is interpolated from pressure_start_
    // to pressure_ using step_progress. An assigned amplitude scales the
    // target pressure at time instead; ignore_amplitude uses the nominal target.
    //
    // The traction at each integration point is t(x) = -p * n(x), where n(x)
    // is the unit normal evaluated from the supplied nodal geometry. The
    // surface integrator handles interpolation and the physical area measure.
    //
    // Only rhs is modified; equations, system_dof_ids and matrix are unchanged.
    void apply(
        model::ModelData&      model_data,
        model::Field&          rhs,
        constraint::Equations& equations,
        const SystemDofIds&    system_dof_ids,
        TripletList&           matrix,
        Precision              time,
        bool                   ignore_amplitude = false,
        Precision              step_progress   = Precision(1)
    ) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
