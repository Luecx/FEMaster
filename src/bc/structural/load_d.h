/**
 * @file load_d.h
 * @brief Defines distributed structural surface-traction loads.
 *
 * `DLoad` prescribes a traction vector over a region of finite-element surfaces.
 * Each surface performs consistent quadrature using its interpolation and
 * physical surface Jacobian so the continuous traction is converted to
 * equivalent nodal generalized forces.
 *
 * @see Condition
 * @see model::SurfaceInterface
 * @see cos::CoordinateSystem
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "../condition.h"
#include "../../constraints/types/equation.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"
#include "../../cos/coordinate_system.h"
#include "../../data/region.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Applies a distributed traction vector to a surface region.
 *
 * For one surface element, the consistent nodal force is
 *
 *     f_e = integral_Gamma_e N^T t dGamma,
 *
 * where `N` is the surface interpolation and `t` is the prescribed traction in
 * global coordinates. If the stored vector is defined in a local basis, the
 * integrand is evaluated as
 *
 *     t(x) = a(t) A(x) t_local,
 *
 * with amplitude `a(t)` and local-to-global basis `A(x)`.
 *
 * `NaN` entries in `values_` denote omitted traction components and are replaced
 * by zero before integration.
 */
struct DLoad : Condition {
    // Types
    using Ptr = std::shared_ptr<DLoad>;

    // Target and initial traction components [tx, ty, tz] per unit surface area.
    // NaN in values_ omits a component, while NaN in values_start_ indicates
    // an unspecified start value. apply() interpolates the prescribed values
    // before integrating the resulting traction over the selected surfaces.
    Vec3 values_       = {NAN, NAN, NAN};
    Vec3 values_start_ = {NAN, NAN, NAN};

    // Target finite-element surfaces over which the traction is integrated
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Optional local basis in which vector components are prescribed
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Construction
    DLoad() = default;
    ~DLoad() override = default;

    // Integrate the prescribed traction over every surface in region_ and
    // accumulate its consistent equivalent nodal forces in the global RHS.
    //
    // Without an amplitude, values_start_ and values_ are interpolated using
    // step_progress. With an amplitude, the target traction is scaled by its
    // value at time instead. ignore_amplitude uses the nominal target traction
    // without scaling or interpolation.
    //
    // NaN target components are omitted. If orientation_ is assigned, the
    // local traction vector is transformed at each integration point because
    // the coordinate-system axes may vary across the surface.
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
