/**
 * @file heat_flux.h
 * @brief Defines prescribed thermal surface heat flux.
 *
 * `HeatFlux` is the thermal Neumann condition. A scalar heat-flow density is
 * integrated over a compiled surface region and distributed consistently to the
 * connected scalar temperature DOFs through the surface shape functions.
 *
 * The condition contributes only to the thermal right-hand side. It does not
 * prescribe temperature and does not modify the conductivity operator.
 *
 * @see Neumann
 * @see ThermalCondition
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include "neumann.h"
#include "../thermal.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Prescribes a uniform scalar heat flux on a surface region.
 *
 * For each selected surface `Gamma_e`, the consistent nodal heat-flow vector is
 *
 *     q_e = integral_Gamma_e N^T q_bar dGamma,
 *
 * where `q_bar` is the amplitude-scaled prescribed heat flux. Positive values
 * follow the thermal RHS sign convention and therefore represent heat input into
 * the assembled balance.
 */
struct HeatFlux : Neumann, ThermalCondition {
    // Types
    using Ptr = std::shared_ptr<HeatFlux>;

    // Compiled surface region receiving the prescribed heat flux
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Nominal heat-flow density per unit reference surface area
    Precision heat_flux_ = Precision(0);

    // Integrate the consistent scalar nodal heat-flow contribution
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
