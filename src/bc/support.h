/******************************************************************************
 * @file support.h
 * @brief Defines the Support class for managing boundary conditions in an FEM model.
 *
 * @details The Support class provides a unified interface for applying constraints to nodes,
 *          elements, or surfaces in an FEM model. Each instance holds a region (node, element,
 *          or surface), a vector of constraint values, and an optional local coordinate system
 *          for transformations. The class ensures proper application of boundary conditions,
 *          with transformations applied where necessary.
 *
 * @date Created on 28.11.2024
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include "../constraints/equation.h"
#include "../core/core.h"
#include "../cos/coordinate_system.h"
#include "../data/region.h"

namespace fem {

/**
 * @struct Support
 * @brief Represents a boundary condition applied to nodes, elements, or surfaces.
 *
 * @details The Support class applies constraints to the specified region (node, element, or surface).
 *          Constraints can be defined in the global or a local coordinate system. The class supports
 *          transformation of constraint values when a local coordinate system is provided.
 */
struct Support {
    using NodeRegionPtr    = model::NodeRegion::Ptr;
    using ElementRegionPtr = model::ElementRegion::Ptr;
    using SurfaceRegionPtr = model::SurfaceRegion::Ptr;

    /// Default constructor
    Support() = default;

    /**
     * @brief Constructs a Support for a node region.
     * @param node_region Region of nodes to apply constraints.
     * @param values Vector of constraint values for the degrees of freedom.
     * @param coordinate_system Optional local coordinate system for transformations.
     */
    Support(NodeRegionPtr node_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Constructs a Support for an element region.
     * @param element_region Region of elements to apply constraints.
     * @param values Vector of constraint values for the degrees of freedom.
     * @param coordinate_system Optional local coordinate system for transformations.
     */
    Support(ElementRegionPtr element_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Constructs a Support for a surface region.
     * @param surface_region Region of surfaces to apply constraints.
     * @param values Vector of constraint values for the degrees of freedom.
     * @param coordinate_system Optional local coordinate system for transformations.
     */
    Support(SurfaceRegionPtr surface_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Applies the support constraints to the model data.
     * @param model_data The FEM model data.
     * @param bc The boundary condition data structure.
     * @param equations The constraint equations structure.
     */
    void apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations);

private:
    NodeRegionPtr node_region = nullptr;       ///< Pointer to the node region, if applicable.
    ElementRegionPtr element_region = nullptr; ///< Pointer to the element region, if applicable.
    SurfaceRegionPtr surface_region = nullptr; ///< Pointer to the surface region, if applicable.
    Vec6 values{NAN, NAN, NAN, NAN, NAN, NAN}; ///< Constraint values for the degrees of freedom.
    cos::CoordinateSystem::Ptr coordinate_system = nullptr; ///< Local coordinate system, if applicable.

    void apply_to_node(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations, ID node_id);
};

} // namespace fem
