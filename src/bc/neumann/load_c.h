/**
 * @file load_c.h
 * @brief Declares concentrated nodal forces and moments.
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

struct CLoad : StructuralNeumann {
    using Ptr = std::shared_ptr<CLoad>;

    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};
    SPtr<model::NodeRegion> region_ = nullptr;

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
