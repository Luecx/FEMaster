/**
 * @file load_p.h
 * @brief Defines scalar pressure loads acting normal to finite-element surfaces.
 *
 * `PLoad` converts a scalar pressure into a traction aligned with the geometric
 * surface normal and delegates consistent finite-element integration, including
 * the analytic follower-load tangent, to the selected surfaces.
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
 *
 * A dead pressure is integrated over the reference configuration and does not
 * depend on the displacement. A follower pressure is integrated over the
 * current configuration, so both the normal and the area measure follow the
 * deformation, and contributes the analytic load tangent df/du.
 */
struct PLoad : Condition {
    // Types
    using Ptr = std::shared_ptr<PLoad>;

    // Nominal pressure per unit surface area, positive against the normal
    Precision pressure_ = NAN;

    // Target finite-element surfaces over which the pressure is integrated
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Construction
    PLoad() = default;
    ~PLoad() override = default;

    // Integrate the pressure traction over every surface in region_ and
    // accumulate the consistent equivalent nodal forces in load. Follower
    // pressures use the current positions and add df/du to load_tangent using
    // the equation indices of system_dof_ids; dead pressures use the reference
    // positions. Constraints and the system matrix remain unchanged.

    void apply(
        model::ModelData&      model_data,
        model::Field&          load,
        TripletList*           load_tangent,
        constraint::Equations& constraints,
        const SystemDofIds&    system_dof_ids,
        TripletList&           system_matrix,
        Precision              scale
    ) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
