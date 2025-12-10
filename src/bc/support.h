/**
 * @file support.h
 * @brief Declares structural supports that constrain degrees of freedom.
 *
 * Supports impose kinematic constraints on selected model entities, producing
 * constraint equations that are later assembled into the solver system.
 *
 * @see src/bc/support.cpp
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../constraints/equation.h"
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
 * @struct Support
 * @brief Represents a kinematic boundary condition over model regions.
 *
 * A support may target node, element, or surface regions. If a local coordinate
 * system is provided, constraint directions are transformed accordingly before
 * being converted into algebraic equations.
 */
struct Support : public fem::Printable {
    using NodeRegionPtr = model::NodeRegion::Ptr;       ///< Alias for node-region pointer.
    using ElementRegionPtr = model::ElementRegion::Ptr; ///< Alias for element-region pointer.
    using SurfaceRegionPtr = model::SurfaceRegion::Ptr; ///< Alias for surface-region pointer.

    /**
     * @brief Default constructor.
     */
    Support() = default;

    /**
     * @brief Creates a support acting on a node region.
     *
     * @param node_region Region of nodes receiving the constraint.
     * @param values Generalized displacement/rotation specification.
     * @param coordinate_system Optional local coordinate system.
     */
    Support(NodeRegionPtr node_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Creates a support acting on an element region.
     *
     * @param element_region Region of elements whose nodes are constrained.
     * @param values Generalized displacement/rotation specification.
     * @param coordinate_system Optional local coordinate system.
     */
    Support(ElementRegionPtr element_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Creates a support acting on a surface region.
     *
     * @param surface_region Region of surfaces whose nodes are constrained.
     * @param values Generalized displacement/rotation specification.
     * @param coordinate_system Optional local coordinate system.
     */
    Support(SurfaceRegionPtr surface_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Applies the support and generates constraint equations.
     *
     * @param model_data FEM model data with geometry and topology.
     * @param equations Container receiving the generated constraint equations.
     */
    void apply(model::ModelData& model_data, constraint::Equations& equations);

    /**
     * @brief One-line description with target and DOF spec.
     */
    std::string str() const override;

private:
    NodeRegionPtr node_region = nullptr;          ///< Targeted node region.
    ElementRegionPtr element_region = nullptr;    ///< Targeted element region.
    SurfaceRegionPtr surface_region = nullptr;    ///< Targeted surface region.
    Vec6 values{NAN, NAN, NAN, NAN, NAN, NAN};    ///< Constraint specification per DOF.
    cos::CoordinateSystem::Ptr coordinate_system = nullptr; ///< Optional local coordinate frame.

    /**
     * @brief Applies the support to a single node and generates equations.
     *
     * @param model_data FEM model data with geometry and topology.
     * @param equations Container receiving the generated constraint equations.
     * @param node_id Identifier of the node to constrain.
     */
    void apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id);
};

} // namespace bc
} // namespace fem
