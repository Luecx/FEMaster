/**
 * @file point_mass.h
 * @brief Point-mass feature contributing diagonal K and M on a node set.
 *
 * Point masses attach translational mass, rotary inertia and optional spring
 * constants directly to existing nodal DOFs.
 *
 * @see src/feature/point_mass.cpp
 * @see src/feature/feature.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "feature.h"

#include "../data/region.h"

namespace fem {
namespace feature {
/**
 * @struct PointMass
 * @brief Concentrated nodal mass, rotary inertia and spring feature.
 *
 * The feature activates every component carrying non-zero mass, inertia or
 * spring data and assembles the corresponding diagonal operator entries.
 */
struct PointMass : Feature {
    using Ptr = std::shared_ptr<PointMass>; ///< Shared pointer alias for point-mass features.

    model::NodeRegion::Ptr region_                  = nullptr;      ///< Target nodes.
    Precision              mass_                    = 0;            ///< Translational mass used for x, y and z.
    Vec3                   rotary_inertia_          = Vec3::Zero(); ///< Rotary inertia on rx, ry and rz.
    Vec3                   spring_constants_        = Vec3::Zero(); ///< Translational spring constants.
    Vec3                   rotary_spring_constants_ = Vec3::Zero(); ///< Rotational spring constants.

    // DOF activation and diagonal operator assembly
    void activate_dofs(SystemDofs& mask) const override;

    /**
     * @brief Adds translational and rotational spring constants to stiffness.
     *
     * @param indices Active global DOF ids. Negative entries are skipped.
     * @param out Triplet list receiving diagonal stiffness entries.
     */
    void assemble_stiffness(const SystemDofIds& indices, TripletList& out) const override;

    /**
     * @brief Adds translational mass and rotary inertia to mass.
     *
     * @param indices Active global DOF ids. Negative entries are skipped.
     * @param out Triplet list receiving diagonal mass entries.
     */
    void assemble_mass(const SystemDofIds& indices, TripletList& out) const override;
};
} // namespace feature
} // namespace fem
