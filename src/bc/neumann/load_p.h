/**
 * @file load_p.h
 * @brief Declares scalar pressure loading on surfaces.
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

struct PLoad : StructuralNeumann {
    using Ptr = std::shared_ptr<PLoad>;

    Precision pressure_ = NAN;
    SPtr<model::SurfaceRegion> region_ = nullptr;

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
