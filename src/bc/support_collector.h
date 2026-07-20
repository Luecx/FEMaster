/**
 * @file support_collector.h
 * @brief Declares the named collection of structural support definitions.
 *
 * `SupportCollector` stores supports by value and offers overloads for node,
 * element and surface targets. When requested, it traverses every stored
 * support and combines their generated algebraic constraints into one equation
 * collection for subsequent solver preprocessing.
 *
 * @see SupportCollector
 * @see Support
 * @see support_collector.cpp
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
 * @brief Stores support definitions and assembles their constraint equations.
 *
 * Supports are held directly in the inherited collection rather than through
 * polymorphic pointers because all target variants share the same concrete
 * `Support` type. The overloaded insertion functions select the corresponding
 * region constructor, while `get_equations()` performs the complete expansion
 * into node-level solver constraints.
 */
struct SupportCollector : model::Collection<Support> {
    // Shared ownership type used by model and load-case registries.
    using Ptr = std::shared_ptr<SupportCollector>;

    // Construct a named support collection with the ordering and uniqueness
    // behavior configured by the generic base constructor.
    explicit SupportCollector(const std::string& name);

    // Destroy the value-stored support entries together with the collection.
    ~SupportCollector() = default;

    // Apply all stored supports to a fresh equation container and return the
    // combined constraints by value.
    constraint::Equations get_equations(model::ModelData& model_data);

    // Provide read-only access to the support definitions in their insertion
    // order.
    const std::vector<Support>& entries() const { return this->_data; }

    // Add a support acting directly on all nodes of `region`.
    void add_supp(model::NodeRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    // Add a support that constrains the nodes connected to every element in
    // `region`.
    void add_supp(model::ElementRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);

    // Add a support that constrains the nodes connected to every surface in
    // `region`.
    void add_supp(model::SurfaceRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system = nullptr);
};
} // namespace bc
} // namespace fem
