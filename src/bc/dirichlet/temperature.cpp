/**
 * @file temperature.cpp
 * @brief Implements prescribed-temperature Dirichlet conditions.
 *
 * A temperature condition may target nodes directly or obtain its nodes from
 * compiled element and surface connectivity. Resolved node identifiers are
 * deduplicated in deterministic traversal order before one scalar thermal
 * constraint equation is emitted per node.
 *
 * The emitted equation uses component zero only. That component is interpreted
 * as the scalar temperature unknown by the dedicated thermal DOF mapping; it is
 * not a structural x-translation despite sharing the same local component index
 * in the generic equation representation.
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
 * only once. The first occurrence of a node determines its position in the
 * generated equation sequence, so traversal remains deterministic.
 *
 * Every equation contains one unit coefficient for thermal component zero and
 * the stored absolute temperature as right-hand side. No structural DOF mapping
 * is consulted here; interpretation of component zero belongs to the thermal
 * constraint collector and its scalar thermal index matrix.
 *
 * @param model_data Compiled model topology and global nodal domain.
 * @param equations Thermal constraint collection receiving the equations.
 */
void Temperature::apply(model::ModelData& model_data, constraint::Equations& equations) {
    // Validate that the semantic definition is unambiguous before any topology
    // expansion. A temperature condition must address exactly one region type.
    const int active_regions = static_cast<int>(node_region_ != nullptr)
                             + static_cast<int>(surface_region_ != nullptr)
                             + static_cast<int>(element_region_ != nullptr);
    logging::error(active_regions == 1,
        "TEMPERATURE: exactly one node, surface or element region must be configured");

    // Reject NaN and infinite prescribed values before they enter the affine
    // constraint right-hand side.
    logging::error(std::isfinite(temperature_),
        "TEMPERATURE: prescribed temperature must be finite");

    // The compiled nodal position field provides the authoritative global node
    // count used to validate all direct and topology-expanded identifiers.
    logging::error(model_data.positions != nullptr,
        "TEMPERATURE: model positions are not initialized");

    // Preserve target traversal order while suppressing nodes shared by several
    // selected elements or surfaces. A dense boolean mask is sufficient because
    // compiled node identifiers are contiguous.
    const Index       node_count = model_data.positions->rows;
    std::vector<bool> selected(node_count, false);
    std::vector<ID>   node_ids{};

    // Centralize node-range validation and deterministic first-occurrence
    // deduplication for all three supported target-region kinds.
    auto add_node = [&](ID node_id) {
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < node_count,
            "TEMPERATURE: node ", node_id, " is outside the compiled node domain");

        if (!selected[static_cast<Index>(node_id)]) {
            selected[static_cast<Index>(node_id)] = true;
            node_ids.push_back(node_id);
        }
    };

    // A direct node region already contains the exact compiled entities carrying
    // the prescribed temperature and requires no topology lookup.
    if (node_region_) {
        for (ID node_id : *node_region_) {
            add_node(node_id);
        }
    }

    // Expand each selected element through its compiled connectivity. Shared
    // corner, edge or higher-order nodes are accepted only on first occurrence.
    if (element_region_) {
        for (ID element_id : *element_region_) {
            // Validate the dense compiled element identifier before resolving the
            // corresponding element pointer.
            logging::error(element_id >= 0 && static_cast<Index>(element_id) < model_data.elements.size(),
                "TEMPERATURE: element ", element_id, " is outside the compiled element domain");

            const auto& element = model_data.elements[static_cast<Index>(element_id)];
            logging::error(element != nullptr,
                "TEMPERATURE: element ", element_id, " is not initialized");

            // Traverse the common element connectivity and feed every attached
            // node through the shared validation and deduplication path.
            for (ID node_id : *element) {
                add_node(node_id);
            }
        }
    }

    // Expand selected surfaces through their compiled connectivity. This is
    // intentionally independent of the parent element because surface regions
    // already represent the exact boundary topology chosen by the user.
    if (surface_region_) {
        for (ID surface_id : *surface_region_) {
            // Validate the dense compiled surface identifier before resolving the
            // corresponding surface pointer.
            logging::error(surface_id >= 0 && static_cast<Index>(surface_id) < model_data.surfaces.size(),
                "TEMPERATURE: surface ", surface_id, " is outside the compiled surface domain");

            const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
            logging::error(surface != nullptr,
                "TEMPERATURE: surface ", surface_id, " is not initialized");

            // Deduplicate nodes shared by adjacent selected faces while keeping
            // the first face/node traversal order stable.
            for (ID node_id : *surface) {
                add_node(node_id);
            }
        }
    }

    // Reserve the final equation capacity once and emit one scalar affine
    // Dirichlet equation for every unique resolved node.
    equations.reserve(equations.size() + node_ids.size());
    for (ID node_id : node_ids) {
        // Component zero denotes temperature only in the thermal system mapping.
        // The unit coefficient fixes that scalar unknown directly to the absolute
        // value stored by this condition.
        const constraint::EquationEntry entry{node_id, Dim(0), Precision(1)};
        equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry}, temperature_);
    }
}

/**
 * Builds the diagnostic representation of the prescribed temperature.
 *
 * The target is reported using the region type and current region cardinality;
 * no compiled connectivity traversal is required for this model-overview text.
 *
 * @return Human-readable target region and absolute temperature.
 */
std::string Temperature::str() const {
    // Resolve the configured region directly from the semantic definition rather
    // than expanding its compiled nodes merely for diagnostic output.
    const auto target = [&]() -> std::string {
        if (node_region_)
            return "NSET " + node_region_->name + " (" + std::to_string(node_region_->size()) + ")";
        if (surface_region_)
            return "SFSET " + surface_region_->name + " (" + std::to_string(surface_region_->size()) + ")";
        if (element_region_)
            return "ELSET " + element_region_->name + " (" + std::to_string(element_region_->size()) + ")";
        return "(unknown)";
    }();

    // Keep the resulting description compact while retaining both target and
    // prescribed scalar value.
    std::ostringstream os;
    os << "TEMPERATURE: target=" << target << ", value=" << temperature_;
    return os.str();
}

} // namespace fem::bc
