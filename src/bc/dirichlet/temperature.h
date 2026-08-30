/**
 * @file temperature.h
 * @brief Declares prescribed-temperature Dirichlet conditions.
 *
 * A `Temperature` assigns one scalar temperature to a node region or to all
 * nodes connected to an element or surface region. Region expansion and
 * construction of scalar thermal constraint equations are implemented in
 * `temperature.cpp`; the boundary-condition object itself stores only the
 * semantic target and prescribed value.
 *
 * The generated equations use local degree of freedom zero as the scalar
 * temperature unknown and therefore belong to the dedicated thermal system
 * mapping rather than to the six-component structural index matrix. Shared
 * nodes encountered through element or surface topology are constrained only
 * once during application.
 *
 * @see Temperature
 * @see Dirichlet
 * @see constraint::Equation
 * @see model::ThermalElement
 * @see temperature.cpp
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "dirichlet.h"

#include "../../core/types_num.h"
#include "../../data/region.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Prescribes one scalar temperature over a compiled model region.
 *
 * Exactly one target region is selected at construction. Node regions are used
 * directly, whereas element and surface regions are expanded through compiled
 * global topology when the condition is applied. A deterministic selection mask
 * suppresses repeated nodes without changing the traversal order of the target
 * definition.
 *
 * `temperature_` is the absolute prescribed temperature and becomes the right-
 * hand side of every generated scalar equation. Each equation contains exactly
 * one unit coefficient for component zero, which is interpreted as temperature
 * only by the thermal DOF mapping used by the thermal loadcase.
 */
struct Temperature : Dirichlet {
    // Shared ownership type for code that needs the concrete prescribed-
    // temperature interface rather than the generic Dirichlet base type.
    using Ptr = std::shared_ptr<Temperature>;

    // Exactly one target region is active for a temperature condition. Element
    // and surface targets are expanded to their compiled nodal connectivity only
    // when the condition is applied.
    model::NodeRegion::Ptr    node_region_    = nullptr;
    model::SurfaceRegion::Ptr surface_region_ = nullptr;
    model::ElementRegion::Ptr element_region_ = nullptr;

    // Absolute scalar temperature prescribed on every unique resolved node.
    Precision temperature_ = Precision(0);

    // Construction from direct or topology-expanded target regions. Each
    // overload records one target kind and leaves the other region pointers null.
    Temperature(model::NodeRegion::Ptr node_region, Precision temperature)
        : node_region_(std::move(node_region)), temperature_(temperature) {}
    Temperature(model::SurfaceRegion::Ptr surface_region, Precision temperature)
        : surface_region_(std::move(surface_region)), temperature_(temperature) {}
    Temperature(model::ElementRegion::Ptr element_region, Precision temperature)
        : element_region_(std::move(element_region)), temperature_(temperature) {}
    ~Temperature() override = default;

    // Resolve the selected target against compiled topology, deduplicate shared
    // nodes and append one scalar thermal Dirichlet equation per unique node.
    void apply(model::ModelData& model_data, constraint::Equations& equations) override;

    // Return the configured target region and prescribed absolute temperature in
    // a compact model-overview representation.
    std::string str() const override;
};

} // namespace fem::bc
