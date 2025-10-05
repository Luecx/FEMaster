/**
 * @file support_collector.h
 * @brief Declares the collector that manages support boundary conditions.
 *
 * The `SupportCollector` aggregates support definitions over multiple regions
 * and converts them into constraint equations on demand.
 *
 * @see src/bc/support_collector.cpp
 * @see src/bc/support.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../data/collection.h"
#include "../data/region.h"
#include "support.h"

namespace fem {
namespace model {
class ModelData;
}
}

namespace fem {
namespace bc {

/**
 * @struct SupportCollector
 * @brief Stores supports and assembles their constraint equations.
 *
 * Supports are stored by value for cache locality. When queried the collector
 * iterates every entry and generates the corresponding constraint equations.
 */
struct SupportCollector : model::Collection<Support> {
    using Ptr = std::shared_ptr<SupportCollector>; ///< Shared pointer alias for collectors.

    /**
     * @brief Constructs a collector with the specified name.
     *
     * @param name Identifier for the collector within the model context.
     */
    explicit SupportCollector(const std::string& name);

    /**
     * @brief Defaulted virtual destructor for polymorphic cleanup.
     */
    ~SupportCollector() = default;

    /**
     * @brief Assembles constraint equations for all stored supports.
     *
     * @param model_data FEM model data used to resolve regions and coordinates.
     * @return constraint::Equations Collection of constraint equations.
     */
    constraint::Equations get_equations(model::ModelData& model_data);

    /**
     * @brief Adds a node-region support to the collector.
     *
     * @param region Node region receiving the support.
     * @param values Constraint specification per generalized DOF.
     * @param coordinate_system Optional local coordinate frame.
     */
    void add_supp(model::NodeRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Adds an element-region support to the collector.
     *
     * @param region Element region receiving the support.
     * @param values Constraint specification per generalized DOF.
     * @param coordinate_system Optional local coordinate frame.
     */
    void add_supp(model::ElementRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    /**
     * @brief Adds a surface-region support to the collector.
     *
     * @param region Surface region receiving the support.
     * @param values Constraint specification per generalized DOF.
     * @param coordinate_system Optional local coordinate frame.
     */
    void add_supp(model::SurfaceRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);
};

} // namespace bc
} // namespace fem
