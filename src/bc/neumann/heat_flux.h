/**
 * @file heat_flux.h
 * @brief Declares prescribed scalar heat flux on surface regions.
 */

#pragma once

#include "thermal_load.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Prescribes a normal heat flux entering the thermal domain.
 *
 * Positive values add heat to the model and therefore contribute positively to
 * the steady-state thermal right-hand side.
 */
struct HeatFlux : ThermalLoad {
    using Ptr = std::shared_ptr<HeatFlux>;

    model::SurfaceRegion::Ptr region_ = nullptr;
    Precision                 flux_   = Precision(0);

    void apply(model::ModelData&   model_data,
               const SystemDofIds& thermal_dof_ids,
               TripletList&        matrix_terms,
               DynamicVector&      heat_source) const override;

    std::string str() const override;
};

} // namespace fem::bc
