/******************************************************************************
 * @file support_collector.h
 * @brief Defines the SupportCollector class for managing collections of supports in an FEM model.
 *
 * @details The SupportCollector class extends the Collection class and provides an interface
 *          to manage and apply multiple Support instances to a FEM model. Supports can be
 *          added for node, element, or surface regions with optional local coordinate systems.
 *
 * @date Created on 27.11.2024
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include "../data/collection.h"
#include "../data/region.h"
#include "support.h"

namespace fem {

/**
 * @struct SupportCollector
 * @brief Manages a collection of Support objects for applying boundary conditions in an FEM model.
 *
 * @details The SupportCollector allows managing multiple supports for various regions (nodes,
 *          elements, and surfaces) and provides functionality to apply all collected supports
 *          to the FEM model. Supports can be added with specific values and optional local
 *          coordinate systems.
 */
struct SupportCollector : model::Collection<Support> {
    using Ptr = std::shared_ptr<SupportCollector>;

    /**
     * @brief Constructor for the SupportCollector.
     * @param p_name Name of the support collector.
     */
    explicit SupportCollector(const std::string& p_name);

    /// Destructor
    ~SupportCollector() = default;

    /**
     * @brief Applies all supports in the collector to the FEM model.
     * @param model_data The FEM model data.
     * @param bc The boundary condition data structure.
     * @param equations The constraint equations structure.
     */
    void apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations);

    /**
     * @brief Adds a Support for a node region.
     * @param region Pointer to the node region.
     * @param values Constraint values for the degrees of freedom.
     * @param coordinate_system Optional local coordinate system for transformations.
     */
    void add_supp(model::NodeRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Adds a Support for an element region.
     * @param region Pointer to the element region.
     * @param values Constraint values for the degrees of freedom.
     * @param coordinate_system Optional local coordinate system for transformations.
     */
    void add_supp(model::ElementRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Adds a Support for a surface region.
     * @param region Pointer to the surface region.
     * @param values Constraint values for the degrees of freedom.
     * @param coordinate_system Optional local coordinate system for transformations.
     */
    void add_supp(model::SurfaceRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);
};

} // namespace fem

