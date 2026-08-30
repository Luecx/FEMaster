/**
 * @file load_inertial.h
 * @brief Declares equivalent rigid-body inertia loading.
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

struct InertialLoad : StructuralNeumann {
    using Ptr = std::shared_ptr<InertialLoad>;

    Vec3 center_     = {0, 0, 0};
    Vec3 center_acc_ = {0, 0, 0};
    Vec3 omega_      = {0, 0, 0};
    Vec3 alpha_      = {0, 0, 0};

    SPtr<model::ElementRegion> region_ = nullptr;
    bool consider_point_masses_ = false;

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
