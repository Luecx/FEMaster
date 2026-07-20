/**
 * @file support.h
 * @brief Declares prescribed structural degrees of freedom over model regions.
 *
 * A `Support` translates displacement and rotation prescriptions into algebraic
 * constraint equations. It can target nodes directly or expand element and
 * surface regions to their connected nodes. Prescriptions may be expressed in
 * the global basis or in a position-dependent coordinate system, in which case
 * one local constraint becomes a linear combination of global degrees of
 * freedom. Region traversal and equation construction are implemented in
 * `support.cpp`.
 *
 * @see Support
 * @see constraint::Equation
 * @see support.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../constraints/types/equation.h"
#include "../core/core.h"
#include "../core/printable.h"
#include "../cos/coordinate_system.h"
#include "../data/region.h"

namespace fem {
namespace model {
class ModelData;
}
}

namespace fem {
namespace bc {

/**
 * @brief Converts regional kinematic prescriptions into solver equations.
 *
 * Exactly one of the node, element or surface region pointers identifies the
 * support target. The six entries of `values_` correspond to three translations
 * followed by three rotations; `NaN` marks a free degree of freedom. In a local
 * coordinate system, the selected local axis is expanded into its three global
 * directional coefficients before the equation is appended.
 */
struct Support : public fem::Printable {
    // Pointer aliases keep the three supported region variants concise in the
    // public constructors and internal state.
    using NodeRegionPtr    = model::NodeRegion::Ptr;
    using ElementRegionPtr = model::ElementRegion::Ptr;
    using SurfaceRegionPtr = model::SurfaceRegion::Ptr;

    // Construct an empty support. This is primarily useful for containers and
    // serializers that assign the target and values later.
    Support() = default;

    // Construct a support that acts directly on every node in `node_region`.
    Support(NodeRegionPtr node_region,
            const Vec6& values,
            cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    // Construct a support that expands each selected element to its connected
    // nodes before generating equations.
    Support(ElementRegionPtr element_region,
            const Vec6& values,
            cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    // Construct a support that expands each selected surface to its connected
    // nodes before generating equations.
    Support(SurfaceRegionPtr surface_region,
            const Vec6& values,
            cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    // Traverse the active target region and append all resulting constraint
    // equations to `equations`. Repeated nodes are processed as encountered by
    // the region topology.
    void apply(model::ModelData& model_data, constraint::Equations& equations);

    // Return a compact representation of the active target, prescribed degrees
    // of freedom and optional coordinate system.
    std::string str() const override;

private:
    // Optional direct node target. When set, no topology expansion is required.
    NodeRegionPtr node_region_ = nullptr;

    // Optional element target. Every node connected to each valid selected
    // element receives the support prescription.
    ElementRegionPtr element_region_ = nullptr;

    // Optional surface target. Every node connected to each valid selected
    // surface receives the support prescription.
    SurfaceRegionPtr surface_region_ = nullptr;

    // Prescribed generalized displacements ordered as
    // `[Ux, Uy, Uz, Rx, Ry, Rz]`. `NaN` leaves the corresponding degree of
    // freedom unconstrained.
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Optional basis in which `values_` is interpreted. A null pointer selects
    // direct constraints in the global coordinate system.
    cos::CoordinateSystem::Ptr coordinate_system_ = nullptr;

    // Generate every active translational or rotational constraint for one
    // node. Global prescriptions produce a single-entry equation; local
    // prescriptions produce three global coefficients for the selected axis.
    void apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id);
};
} // namespace bc
} // namespace fem
