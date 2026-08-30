/**
 * @file thermal_load.h
 * @brief Declares the common interface for steady-state thermal Neumann data.
 */

#pragma once

#include "neumann.h"
#include "../../core/types_eig.h"

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Natural thermal boundary condition assembled into K_T and f_T.
 */
struct ThermalLoad : Neumann {
    using Ptr = std::shared_ptr<ThermalLoad>;

    ~ThermalLoad() override = default;

    virtual void apply(model::ModelData&  model_data,
                       const SystemDofIds& thermal_dof_ids,
                       TripletList&        matrix_terms,
                       DynamicVector&      heat_source) const = 0;
};

} // namespace fem::bc
