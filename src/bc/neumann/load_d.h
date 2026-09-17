/**
 * @file load_d.h
 * @brief Defines distributed surface-traction Neumann conditions.
 *
 * `DLoad` prescribes a traction vector over a region of finite-element surfaces.
 * Each surface performs consistent quadrature using its interpolation and
 * physical surface Jacobian so the continuous traction is converted to
 * equivalent nodal generalized forces.
 *
 * @see Neumann
 * @see model::SurfaceInterface
 * @see cos::CoordinateSystem
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
struct DLoad : Neumann {
    // Types
    using Ptr = std::shared_ptr<DLoad>;

    // Nominal traction components per physical surface area. NaN marks an
    // omitted component.
    Vec3 values_ = {NAN, NAN, NAN};

    // Target finite-element surfaces over which the traction is integrated
    SPtr<model::SurfaceRegion> region_ = nullptr;

    // Construction
    DLoad() = default;
    ~DLoad() override = default;

    // Integrate N^T t over every selected physical surface and accumulate the
    // equivalent nodal forces in the supplied RHS field.
    void apply(model::ModelData& model_data, model::Field& rhs,
               Precision time, bool ignore_amplitude = false) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
