/**
 * @file support.h
 * @brief Declares prescribed structural degrees of freedom over model regions.
 *
 * A `Support` is a structural Dirichlet condition that translates displacement
 * and rotation prescriptions into algebraic constraint equations. It can target
 * nodes directly or expand element and surface regions to their connected nodes.
 * Prescriptions may be expressed in the global basis or in a position-dependent
 * coordinate system.
 *
 * @see Support
 * @see Dirichlet
 * @see constraint::Equation
 * @see support.cpp
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "dirichlet.h"

#include "../../core/core.h"
#include "../../cos/coordinate_system.h"
#include "../../data/region.h"

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Converts regional kinematic prescriptions into solver equations.
 *
 * Exactly one of the node, element or surface region pointers identifies the
 * support target. The six entries of `values_` correspond to three translations
 * followed by three rotations; `NaN` marks a free degree of freedom. In a local
 * coordinate system, the selected local axis is expanded into its three global
 * directional coefficients before the equation is appended.
 */
struct Support : Dirichlet {
    using NodeRegionPtr    = model::NodeRegion::Ptr;
    using ElementRegionPtr = model::ElementRegion::Ptr;
    using SurfaceRegionPtr = model::SurfaceRegion::Ptr;

private:
    // Exactly one target region is active for a support definition.
    NodeRegionPtr    node_region_    = nullptr;
    ElementRegionPtr element_region_ = nullptr;
    SurfaceRegionPtr surface_region_ = nullptr;

    // Prescribed generalized displacements ordered as
    // [Ux, Uy, Uz, Rx, Ry, Rz]. NaN leaves the corresponding DOF free.
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Optional basis in which the prescribed components are interpreted.
    cos::CoordinateSystem::Ptr coordinate_system_ = nullptr;

public:
    // Construction from direct or topology-expanded structural targets.
    Support() = default;
    Support(NodeRegionPtr node_region,
            const Vec6& values,
            cos::CoordinateSystem::Ptr coordinate_system = nullptr);
    Support(ElementRegionPtr element_region,
            const Vec6& values,
            cos::CoordinateSystem::Ptr coordinate_system = nullptr);
    Support(SurfaceRegionPtr surface_region,
            const Vec6& values,
            cos::CoordinateSystem::Ptr coordinate_system = nullptr);
    ~Support() override = default;

    // Expand the selected target and append structural Dirichlet equations.
    void apply(model::ModelData& model_data, constraint::Equations& equations) override;

    // Return the target, prescribed DOFs and optional coordinate system.
    std::string str() const override;

private:
    // Generate all active translational or rotational equations for one node.
    void apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id);
};

} // namespace fem::bc
