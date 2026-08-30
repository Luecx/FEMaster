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
 * The ambient term is assembled directly into the scalar nodal RHS. The
 * temperature-dependent term is emitted as symbolic Robin rows and remains in
 * nodal coordinates until Model assembles the thermal system.
 */
struct Convection : Robin {
    using Ptr = std::shared_ptr<Convection>;

    model::SurfaceRegion::Ptr region_ = nullptr;
    Precision film_coefficient_ = Precision(0);
    Precision ambient_temperature_ = Precision(0);

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               RobinEquations&   equations,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
