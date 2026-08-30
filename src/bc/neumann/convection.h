/**
 * @file convection.h
 * @brief Declares thermal convection boundary conditions.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

/**
 * @brief Applies a linear convection boundary to a surface region.
 *
 * The ambient contribution `h * T_inf` is assembled through the ordinary
 * Neumann `apply()` path into component zero of the three-component thermal load
 * field. The temperature-dependent Robin contribution is assembled separately
 * into the thermal conductivity matrix.
 */
struct Convection : Neumann {
    using Ptr = std::shared_ptr<Convection>;

    model::SurfaceRegion::Ptr region_ = nullptr;
    Precision film_coefficient_ = Precision(0);
    Precision ambient_temperature_ = Precision(0);

    Convection() = default;
    ~Convection() override = default;

    void apply(model::ModelData& model_data,
               model::Field&     bc,
               Precision         time,
               bool              ignore_amplitude = false) override;

    void apply_conductivity(model::ModelData& model_data,
                            const SystemDofIds& thermal_dof_ids,
                            TripletList&         triplets,
                            Precision            time,
                            bool                 ignore_amplitude = false);

    std::string str() const override;
};

} // namespace fem::bc
