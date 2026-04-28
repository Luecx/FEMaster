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
    using Ptr = std::shared_ptr<InertialLoad>; ///< Shared pointer alias for inertial loads.

    Vec3                       center_                = {0, 0, 0}; ///< Reference point for rotational terms.
    Vec3                       center_acc_            = {0, 0, 0}; ///< Translational acceleration at the reference point.
    Vec3                       omega_                 = {0, 0, 0}; ///< Angular velocity vector.
    Vec3                       alpha_                 = {0, 0, 0}; ///< Angular acceleration vector.
    SPtr<model::ElementRegion> region_                = nullptr;   ///< Target element region.
    bool                       consider_point_masses_ = false;     ///< Include point-mass features in inertial force assembly.

    /**
     * @brief Default constructor for delayed field assignment by collectors.
     */
    InertialLoad() = default;

    /**
     * @brief Defaulted virtual destructor.
     */
    ~InertialLoad() override = default;

    /**
     * @brief Integrates rigid-body inertia over target structural elements.
     *
     * @param model_data Model data used for elements, nodes and optional point masses.
     * @param bc Boundary-condition field receiving equivalent nodal forces.
     * @param time Current analysis time. Inertial loads currently ignore it.
     */
    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;

    /**
     * @brief Returns a compact description of the inertial load.
     *
     * @return std::string Target region and rigid-body acceleration parameters.
     */
    std::string str() const override;
};
} // namespace bc
} // namespace fem
