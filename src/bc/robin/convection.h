/**
 * @file convection.h
 * @brief Declares linear thermal convection as a Robin boundary condition.
 *
 * `Convection` couples the surface temperature to a prescribed ambient
 * temperature through a non-negative film coefficient. In weak form the
 * ambient term contributes to the scalar thermal load vector, while the part
 * proportional to the unknown surface temperature contributes a symmetric
 * boundary operator.
 *
 * The condition is integrated on compiled surface regions in the reference
 * configuration. Ordinary scalar surface integration is used for the ambient
 * source, whereas a dedicated shape-product integration routine evaluates the
 * `N N^T` operator term with sufficient quadrature order.
 *
 * @see Convection
 * @see Robin
 * @see model::SurfaceInterface
 * @see convection.cpp
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "robin.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Applies linear convection `q = h (T_inf - T)` on a surface region.
 *
 * The thermal weak form separates naturally into
 *
 * \f[
 *     f_h^e = \int_{\Gamma_e} h\,T_\infty N\,\mathrm{d}\Gamma
 * \f]
 *
 * and
 *
 * \f[
 *     K_h^e = \int_{\Gamma_e} h\,N N^T\,\mathrm{d}\Gamma.
 * \f]
 *
 * `f_h^e` is accumulated into the one-component nodal thermal RHS. `K_h^e` is
 * mapped through the supplied scalar temperature DOF matrix and appended to the
 * global sparse triplet list. The resulting global equation therefore contains
 * `K_T + K_h` on the left-hand side.
 *
 * The optional amplitude scales the film coefficient rather than the ambient
 * temperature. A zero effective coefficient makes the condition inactive.
 */
struct Convection : Robin {
    // Shared ownership type for direct references to convection conditions.
    using Ptr = std::shared_ptr<Convection>;

    // Compiled surface region exchanging heat with the ambient environment.
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Nominal non-negative film coefficient h. Optional amplitude scaling is
    // applied to this coefficient during assembly.
    Precision film_coefficient_ = Precision(0);

    // Prescribed ambient temperature T_inf appearing only in the source part of
    // the Robin condition.
    Precision ambient_temperature_ = Precision(0);

    // Assemble the ambient source into `rhs` and the temperature-dependent
    // boundary operator into `lhs` using `system_dof_ids` for nodal mapping.
    void apply(model::ModelData&   model_data,
               model::Field&       rhs,
               const SystemDofIds& system_dof_ids,
               TripletList&         lhs,
               Precision           time,
               bool                ignore_amplitude = false) override;

    // Return the target surface set, nominal film coefficient, ambient
    // temperature and optional amplitude in a compact diagnostic string.
    std::string str() const override;
};

} // namespace fem::bc
