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
 * The class stores only the semantic support definition. Topology expansion and
 * construction of concrete `constraint::Equation` objects occur when the
 * support is applied to compiled model data, so no duplicate node connectivity
 * is retained in the boundary-condition object itself.
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
 *
 * Node regions are already solver-facing targets. Element and surface regions
 * retain their input-level semantics and are expanded through compiled
 * connectivity only during `apply()`. The support does not own or copy the
 * referenced regions or coordinate system beyond shared ownership pointers.
 */
struct Support : Dirichlet {
    // Target-region pointer types retained here for readable constructors and to
    // make the three mutually exclusive support targets explicit in the class.
    using NodeRegionPtr    = model::NodeRegion::Ptr;
    using ElementRegionPtr = model::ElementRegion::Ptr;
    using SurfaceRegionPtr = model::SurfaceRegion::Ptr;

private:
    // Exactly one target region is active for a support definition. Element and
    // surface targets are expanded to their connected nodes when equations are
    // generated against compiled model data.
    NodeRegionPtr    node_region_    = nullptr;
    ElementRegionPtr element_region_ = nullptr;
    SurfaceRegionPtr surface_region_ = nullptr;

    // Prescribed generalized displacements ordered as [Ux, Uy, Uz, Rx, Ry, Rz].
    // NaN is the explicit sentinel for an unconstrained component and therefore
    // produces no equation.
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Optional basis in which the six stored components are interpreted. A null
    // pointer means that translations and rotations already refer to the global
    // model axes.
    cos::CoordinateSystem::Ptr coordinate_system_ = nullptr;

public:
    // Construction from direct or topology-expanded structural targets. The
    // overload itself records the target kind; exactly one region pointer becomes
    // active and no connectivity is expanded at construction time.
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

    // Resolve the selected region through compiled topology and append one
    // structural Dirichlet equation for every finite prescribed generalized
    // component of every encountered node.
    void apply(model::ModelData& model_data, constraint::Equations& equations) override;

    // Return the active target, prescribed generalized DOFs and optional local
    // coordinate system in a compact model-overview representation.
    std::string str() const override;

private:
    // Generate all active translational or rotational equations for one compiled
    // node. Local prescriptions are expanded into three global coefficients using
    // the coordinate-system basis evaluated at that node position.
    void apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id);
};

} // namespace fem::bc
