/**
 * @file convection.h
 * @brief Declares linear thermal convection as a Robin boundary condition.
 */

#pragma once

#include "robin.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Applies q = h (T_inf - T) on a surface region.
 *
 * The ambient term contributes to the scalar thermal RHS. The
 * temperature-dependent term is integrated as
 *
 *     K_h = integral_Gamma h N N^T dGamma
 *
 * and assembled directly into the global thermal operator through the supplied
 * temperature DOF map.
 */
struct Convection : Robin {
    using Ptr = std::shared_ptr<Convection>;

    model::SurfaceRegion::Ptr region_ = nullptr;
    Precision film_coefficient_ = Precision(0);
    Precision ambient_temperature_ = Precision(0);

    void apply(model::ModelData&  model_data,
               model::Field&      rhs,
               const SystemDofIds& system_dof_ids,
               TripletList&        lhs,
               Precision           time,
               bool                ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
