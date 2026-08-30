/**
 * @file heat_flux.h
 * @brief Declares prescribed thermal surface heat flux.
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

struct HeatFlux : ThermalNeumann {
    using Ptr = std::shared_ptr<HeatFlux>;

    model::SurfaceRegion::Ptr region_ = nullptr;
    Precision heat_flux_ = Precision(0);

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
