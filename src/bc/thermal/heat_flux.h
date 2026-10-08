/**
 * @file heat_flux.h
 * @brief Defines prescribed thermal surface heat flux.
 *
 * `HeatFlux` integrates a prescribed scalar surface heat-flow density into
 * equivalent nodal thermal loads using the reference surface geometry.
 * It contributes only to the thermal RHS; no temperatures are prescribed
 * and the conductivity matrix remains unchanged.
 *
 * @see Condition
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include "../condition.h"
#include "../../constraints/types/equation.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"
#include "../../data/region.h"

#include <cmath>

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
struct HeatFlux : Condition {
    // Types
    using Ptr = std::shared_ptr<HeatFlux>;

    // Compiled surfaces on which the prescribed heat flux is integrated
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Target and initial heat-flow densities per unit reference surface area.
    // A NaN initial value indicates an unspecified start, interpreted as zero
    // when interpolating the prescribed heat flux over a step.
    Precision heat_flux_       = Precision(0);
    Precision heat_flux_start_ = NAN;

    // Integrate the prescribed scalar heat flux over every surface in region_
    // and accumulate equivalent nodal heat-flow contributions in rhs.
    //
    // Without an amplitude, the effective flux is interpolated from
    // heat_flux_start_ to heat_flux_ using step_progress. An assigned amplitude
    // instead scales the target flux at time. ignore_amplitude uses the nominal
    // target flux without scaling or interpolation.
    //
    // Each surface integrates q_e = integral(N^T * q dGamma) in the reference
    // configuration, including its shape functions and physical area Jacobian.
    // Positive flux adds heat to the scalar thermal balance.
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
