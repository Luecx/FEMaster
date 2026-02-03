/**
 * @file point_mass.h
 * @brief Point-mass feature contributing diagonal K and M on a node set.
 */

#pragma once

#include "feature.h"
#include "../data/region.h"

namespace fem { namespace feature {

struct PointMass : Feature {
    using Ptr = std::shared_ptr<PointMass>;

    model::NodeRegion::Ptr region;        // target nodes
    Precision              mass = 0;      // translational mass
    Vec3                   rotary_inertia = Vec3::Zero();
    Vec3                   spring_constants = Vec3::Zero();
    Vec3                   rotary_spring_constants = Vec3::Zero();

    void assemble_stiffness(const SystemDofIds& indices, TripletList& out) const override;
    void assemble_mass     (const SystemDofIds& indices, TripletList& out) const override;
};

} } // namespace fem::feature

