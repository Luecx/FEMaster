/**
 * @file thermal_temp.h
 * @brief Declares prescribed-temperature boundary conditions.
 *
 * A `Temperature` assigns one scalar temperature to a node region or to all
 * nodes connected to an element or surface region. Region expansion and
 * construction of scalar thermal constraint equations are implemented in
 * `thermal_temp.cpp`.
 *
 * The generated equations use local degree of freedom zero as the scalar
 * temperature unknown. They belong to a thermal system and require a dedicated
 * thermal DOF mapping rather than the structural six-DOF mapping.
 *
 * @see Temperature
 * @see constraint::Equation
 * @see model::ThermalElement
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#pragma once

#include "bc.h"

#include "../constraints/types/equation.h"
#include "../core/printable.h"
#include "../core/types_num.h"
#include "../data/region.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Prescribes one scalar temperature over a model region.
 *
 * Exactly one target region is selected at construction. Node regions are used
 * directly, whereas element and surface regions are expanded through compiled
 * global topology when the condition is applied. Shared nodes are constrained
 * only once.
 *
 * `temperature_` is the absolute prescribed temperature and becomes the right-
 * hand side of every generated scalar equation. Local degree of freedom zero
 * represents the thermal unknown; the resulting equations are not structural
 * displacement constraints.
 */
struct Temperature : BoundaryCondition, Printable {
    // Shared ownership type for prescribed-temperature definitions
    using Ptr = std::shared_ptr<Temperature>;

    // Optional target regions. Construction initializes exactly one pointer.
    model::NodeRegion::Ptr    node_region_    = nullptr;
    model::SurfaceRegion::Ptr surface_region_ = nullptr;
    model::ElementRegion::Ptr element_region_ = nullptr;

    // Absolute temperature prescribed on every resolved node
    Precision temperature_ = Precision(0);

    // Construction from direct or topology-expanded target regions
    Temperature(model::NodeRegion::Ptr node_region, Precision temperature)
        : node_region_(std::move(node_region)),
          temperature_(temperature) {}

    Temperature(model::SurfaceRegion::Ptr surface_region, Precision temperature)
        : surface_region_(std::move(surface_region)),
          temperature_(temperature) {}

    Temperature(model::ElementRegion::Ptr element_region, Precision temperature)
        : element_region_(std::move(element_region)),
          temperature_(temperature) {}

    ~Temperature() override = default;

    // Resolve the target to unique compiled nodes and append scalar thermal
    // Dirichlet equations with `temperature_` as their right-hand side.
    void apply(model::ModelData& model_data, constraint::Equations& equations);

    // Return the target region and prescribed scalar temperature.
    std::string str() const override;
};

} // namespace fem::bc
