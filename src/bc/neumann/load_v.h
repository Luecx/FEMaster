/**
 * @file load_v.h
 * @brief Declares distributed body-force loading.
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

struct VLoad : StructuralNeumann {
    using Ptr = std::shared_ptr<VLoad>;

    Vec3 values_ = {NAN, NAN, NAN};
    SPtr<model::ElementRegion> region_ = nullptr;

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
