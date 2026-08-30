/**
 * @file heat_flux.h
 * @brief Declares prescribed surface heat-flux boundary conditions.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

/**
 * @brief Applies a uniform prescribed heat flux to a surface region.
 *
 * Positive `heat_flux_` denotes heat entering the model. Thermal Neumann loads
 * reuse the existing three-component nodal load integration and store their
 * scalar contribution exclusively in component zero.
 */
struct HeatFlux : Neumann {
    using Ptr = std::shared_ptr<HeatFlux>;

    model::SurfaceRegion::Ptr region_ = nullptr;
    Precision heat_flux_ = Precision(0);

    HeatFlux() = default;
    ~HeatFlux() override = default;

    void apply(model::ModelData& model_data,
               model::Field&     bc,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
