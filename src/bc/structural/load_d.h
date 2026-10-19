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

    // Nominal generalized components. NaN omits a component.
    Vec3 values_ = {NAN, NAN, NAN};

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
    // Each step evaluates this nominal traction or force density with its
    // initial or target amplitude factor; the load itself remains unchanged.
    //
    // NaN target components are omitted. If orientation_ is assigned, the
    // local traction vector is transformed at each integration point because
    // the coordinate-system axes may vary across the surface.
    //
    // Only rhs is modified; equations, system_dof_ids and matrix are unchanged.
    // Construct an independent physical condition for a new step.
    Condition::Ptr clone() const override { return std::make_shared<DLoad>(*this); }

    void apply(
        model::ModelData&      model_data,
        model::Field&          rhs,
        constraint::Equations& equations,
        const SystemDofIds&    system_dof_ids,
        TripletList&           matrix,
        Precision              time,
        bool                   ignore_amplitude = false,
        Precision              step_progress = Precision(1),
        const LoadInterpolation* interpolation = nullptr
    ) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
