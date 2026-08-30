/**
 * @file load_d.h
 * @brief Declares distributed vector tractions on surfaces.
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

struct DLoad : StructuralNeumann {
    using Ptr = std::shared_ptr<DLoad>;

    Vec3 values_ = {NAN, NAN, NAN};
    SPtr<model::SurfaceRegion> region_ = nullptr;

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
