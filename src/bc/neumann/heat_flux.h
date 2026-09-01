/**
 * @file heat_flux.h
 * @brief Declares prescribed thermal surface heat flux.
 *
 * `HeatFlux` represents a scalar Neumann boundary condition on a compiled
 * surface region. The prescribed heat flux is integrated with the surface shape
 * functions over the reference surface and accumulated as a consistent scalar
 * nodal heat-flow vector. An optional amplitude scales the nominal value at the
 * current analysis time.
 *
 * The class contributes only to the thermal right-hand side. It does not modify
 * the conductivity matrix and therefore remains distinct from convection and
 * other Robin conditions.
 *
 * @see HeatFlux
 * @see ThermalNeumann
 * @see model::SurfaceInterface
 * @see heat_flux.cpp
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Prescribes a uniform scalar heat flux over a surface region.
 *
 * The nominal `heat_flux_` is interpreted as heat flow per unit reference area.
 * For every selected surface, the consistent nodal contribution is
 *
 * \f[
 *     f_q^e = \int_{\Gamma_e} N^T q\,\mathrm{d}\Gamma,
 * \f]
 *
 * where `N` denotes the surface interpolation vector. The surface abstraction
 * owns quadrature, physical area scaling and nodal scattering. This condition
 * only supplies the scalar field value and optional amplitude multiplier.
 */
struct HeatFlux : ThermalNeumann {
    // Shared ownership type for direct references to prescribed heat-flux loads.
    using Ptr = std::shared_ptr<HeatFlux>;

    // Compiled surface region over which the scalar heat flux is integrated.
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Nominal heat flux per unit reference surface area. Positive values follow
    // the thermal RHS sign convention used by the conduction balance.
    Precision heat_flux_ = Precision(0);

    // Evaluate the optional amplitude and integrate the resulting scalar heat
    // flux consistently into the one-component nodal thermal RHS field.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    // Return the target surface set, nominal flux and optional amplitude in a
    // compact diagnostic representation.
    std::string str() const override;
};

} // namespace fem::bc
