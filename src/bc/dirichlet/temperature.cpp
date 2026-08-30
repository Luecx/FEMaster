/**
 * @file temperature.cpp
 * @brief Implements prescribed-temperature Dirichlet conditions.
 *
 * A temperature condition may target nodes directly or obtain its nodes from
 * compiled element and surface connectivity. Resolved node identifiers are
 * deduplicated in deterministic traversal order before one scalar thermal
 * constraint equation is emitted per node.
 *
 * @see Temperature
 * @see constraint::Equation
 * @see model::ThermalElement
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "temperature.h"

#include "../../core/logging.h"
#include "../../model/element/element.h"
#include "../../model/geometry/surface/surface_interface.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <vector>

namespace fem::bc {

/**
 * Resolves the active target region and appends nodal temperature constraints.
 *
 * Exactly one node, element or surface region must be configured. Topological
 * targets are expanded through the compiled model and shared nodes are emitted
 * only once. Every equation contains one unit coefficient for thermal component
 * zero and the stored absolute temperature as right-hand side.
 *
 * @param model_data Compiled model topology and global nodal domain.
 * @param equations Thermal constraint collection receiving the equations.
 */
void Temperature::apply(model::ModelData& model_data, constraint::Equations& equations) {
    // Validate the physical prescription and compiled nodal domain.
    const int active_regions = static_cast<int>(node_region_ != nullptr)
                             + static_cast<int>(surface_region_ != nullptr)
                             + static_cast<int>(element_region_ != nullptr);
    logging::error(active_regions == 1,
        "TEMPERATURE: exactly one node, surface or element region must be configured");
    logging::error(std::isfinite(temperature_),
        "TEMPERATURE: prescribed temperature must be finite");
    logging::error(model_data.positions != nullptr,
        "TEMPERATURE: model positions are not initialized");

    // Preserve target traversal order while suppressing shared nodes.
    const Index       node_count = model_data.positions->rows;
    std::vector<bool> selected(node_count, false);
    std::vector<ID>   node_ids{};

    auto add_node = [&](ID node_id) {
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < node_count,
            "TEMPERATURE: node ", node_id, " is outside the compiled node domain");

        if (!selected[static_cast<Index>(node_id)]) {
            selected[static_cast<Index>(node_id)] = true;
            node_ids.push_back(node_id);
        }
    };

    // Resolve direct node targets without topology traversal.
    if (node_region_) {
        for (ID node_id : *node_region_) {
            add_node(node_id);
        }
    }

    // Expand selected elements to their compiled connectivity.
    if (element_region_) {
        for (ID element_id : *element_region_) {
            logging::error(element_id >= 0 && static_cast<Index>(element_id) < model_data.elements.size(),
                "TEMPERATURE: element ", element_id, " is outside the compiled element domain");

            const auto& element = model_data.elements[static_cast<Index>(element_id)];
            logging::error(element != nullptr,
                "TEMPERATURE: element ", element_id, " is not initialized");

            for (ID node_id : *element) {
                add_node(node_id);
            }
        }
    }

    // Expand selected surfaces and deduplicate nodes shared by adjacent faces.
    if (surface_region_) {
        for (ID surface_id : *surface_region_) {
            logging::error(surface_id >= 0 && static_cast<Index>(surface_id) < model_data.surfaces.size(),
                "TEMPERATURE: surface ", surface_id, " is outside the compiled surface domain");

            const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
            logging::error(surface != nullptr,
                "TEMPERATURE: surface ", surface_id, " is not initialized");

            for (ID node_id : *surface) {
                add_node(node_id);
            }
        }
    }

    // Emit one scalar temperature equation for every unique resolved node.
    equations.reserve(equations.size() + node_ids.size());
    for (ID node_id : node_ids) {
        const constraint::EquationEntry entry{node_id, Dim(0), Precision(1)};
        equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry}, temperature_);
    }
}

/**
 * Builds the diagnostic representation of the prescribed temperature.
 *
 * @return Human-readable target region and absolute temperature.
 */
std::string Temperature::str() const {
    // Resolve the configured region without traversing compiled topology.
    const auto target = [&]() -> std::string {
        if (node_region_)
            return "NSET " + node_region_->name + " (" + std::to_string(node_region_->size()) + ")";
        if (surface_region_)
            return "SFSET " + surface_region_->name + " (" + std::to_string(surface_region_->size()) + ")";
        if (element_region_)
            return "ELSET " + element_region_->name + " (" + std::to_string(element_region_->size()) + ")";
        return "(unknown)";
    }();

    std::ostringstream os;
    os << "TEMPERATURE: target=" << target << ", value=" << temperature_;
    return os.str();
}

} // namespace fem::bc
