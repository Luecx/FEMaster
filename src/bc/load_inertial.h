/**
 * @file load_inertial.h
 * @brief Declares the rigid-body inertial load boundary condition.
 *
 * @see src/bc/load_inertial.cpp
 * @see src/bc/load.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "load.h"

#include "../data/region.h"

namespace fem {
namespace bc {

/**
 * @struct InertialLoad
 * @brief Rigid-body inertial load with translational and rotational terms.
 *
 * The equivalent nodal forces are computed by integrating
 * `f(x) = -rho * a(x)` over the selected element region, where
 * `a(x) = a0 + alpha x r + omega x (omega x r)` and `r = x - center`.
 */
struct InertialLoad : public Load {
    using Ptr = std::shared_ptr<InertialLoad>;

    Vec3 center{0, 0, 0};
    Vec3 center_acc{0, 0, 0};
    Vec3 omega{0, 0, 0};
    Vec3 alpha{0, 0, 0};
    SPtr<model::ElementRegion> region = nullptr; ///< Target element region.
    bool consider_point_masses = false;          ///< Include point-mass features in inertial force assembly.

    InertialLoad() = default;
    ~InertialLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;
    std::string str() const override;
};

} // namespace bc
} // namespace fem
