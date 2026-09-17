/**
 * @file support.h
 * @brief Defines structural Dirichlet conditions on generalized nodal DOFs.
 *
 * A structural support prescribes selected translational and rotational nodal
 * components. Targets may be specified directly as node regions or indirectly
 * through element and surface regions, which are expanded to their connected
 * nodes when the condition is applied.
 *
 * In global coordinates, a prescribed generalized component produces a scalar
 * equation with one unit coefficient. With a local coordinate system, the
 * selected local axis is projected onto the corresponding three global
 * components, producing a row of the structural constraint matrix.
 *
 * @see Dirichlet
 * @see constraint::Equation
 * @see cos::CoordinateSystem
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "dirichlet.h"
#include "../../core/core.h"
#include "../../cos/coordinate_system.h"
#include "../../data/region.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Prescribes structural translations or rotations on a model region.
 *
 * The six entries of `values_` follow the generalized nodal ordering
 *
 *     [Ux, Uy, Uz, Rx, Ry, Rz].
 *
 * A finite entry is prescribed and `NaN` denotes a free component. Exactly one
 * of the node, element or surface region pointers is expected to define the
 * target. Element and surface targets are converted to nodal prescriptions by
 * traversing their compiled connectivity.
 *
 * Without `coordinate_system_`, component `i` generates
 *
 *     u_i = value_i.
 *
 * For an oriented translational component with local unit axis `e_a`, the
 * generated equation is
 *
 *     e_a,x Ux + e_a,y Uy + e_a,z Uz = value_a,
 *
 * and the same projection is applied to `[Rx, Ry, Rz]` for rotational
 * components. The local basis may depend on the nodal position and is therefore
 * evaluated separately for every target node.
 */
struct Support : Dirichlet {
    // Types
    using Ptr              = std::shared_ptr<Support>;
    using NodeRegionPtr    = model::NodeRegion::Ptr;
    using ElementRegionPtr = model::ElementRegion::Ptr;
    using SurfaceRegionPtr = model::SurfaceRegion::Ptr;

private:
    // Exactly one target region defines the entities whose connected nodes are
    // constrained.
    NodeRegionPtr    node_region_    = nullptr;
    ElementRegionPtr element_region_ = nullptr;
    SurfaceRegionPtr surface_region_ = nullptr;

    // Prescribed generalized nodal values ordered as
    // [Ux, Uy, Uz, Rx, Ry, Rz]. NaN marks an unconstrained component.
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Optional local basis in which the six prescribed components are defined
    cos::CoordinateSystem::Ptr coordinate_system_ = nullptr;

public:
    // Construction from node, element or surface targets
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

    // Expand the active region and append one scalar constraint equation for
    // every finite prescribed generalized component at every resolved node.
    void apply(model::ModelData& model_data, constraint::Equations& equations) override;

    // Diagnostics
    std::string str() const override;

private:
    // Convert the six stored prescriptions for one global node into rows of
    // C u = d, including local-to-global projection when an orientation exists.
    void apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id);
};

} // namespace fem::bc
