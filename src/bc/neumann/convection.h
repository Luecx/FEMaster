/**
 * @file convection.h
 * @brief Declares linear film-convection boundary conditions.
 */

#pragma once

#include "thermal_load.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Exchanges heat with an ambient temperature through a film coefficient.
 */
struct Convection : ThermalLoad {
    using Ptr = std::shared_ptr<Convection>;

    model::SurfaceRegion::Ptr region_              = nullptr;
    Precision                 film_coefficient_    = Precision(0);
    Precision                 ambient_temperature_ = Precision(0);

    void apply(model::ModelData&   model_data,
               const SystemDofIds& thermal_dof_ids,
               TripletList&        matrix_terms,
               DynamicVector&      heat_source) const override;

    std::string str() const override;
};

} // namespace fem::bc
