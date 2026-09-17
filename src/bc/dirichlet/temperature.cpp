/**
 * @file temperature.cpp
 * @brief Implements region expansion and scalar temperature constraints.
 *
 * A prescribed temperature is a scalar Dirichlet condition. The semantic target
 * may already be nodal or may require expansion through compiled element or
 * surface connectivity. Shared nodes are deduplicated before one equation
 *
 *     T_i = T_bar
 *
 * is appended per physical node.
 *
 * The implementation does not assemble heat flux, conductivity or convection
 * terms. It only constructs the algebraic primary-variable prescriptions that a
 * thermal analysis applies to its scalar temperature system.
 *
 * @see Temperature
 * @see Dirichlet
 * @see constraint::Equation
 *
 * @author Finn Eggers
 * @date 17.09.2026
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
 * Resolves the semantic target and appends nodal temperature equations.
 *
 * Exactly one target region must be configured. Node regions are consumed
 * directly. Element and surface regions are expanded through the compiled
 * connectivity and all resulting node identifiers are deduplicated in first-
 * encounter order.
 *
 * With scalar thermal DOF zero, every unique node `i` contributes the equation
 *
 *     1 * T_i = T_bar,
 *
 * where `T_bar` is the stored absolute temperature. No modification of the
 * thermal conductivity matrix or load vector occurs here.
 *
 * @param model_data Compiled topology and nodal domain used for target expansion.
 * @param equations Constraint collection receiving one scalar equation per
 *                  unique target node.
 */
void Temperature::apply(model::ModelData& model_data, constraint::Equations& equations) {
    // Validate the semantic definition and the compiled nodal domain before any
    // topology is traversed
    const int active_regions = static_cast<int>(node_region_    != nullptr)
                             + static_cast<int>(surface_region_ != nullptr)
                             + static_cast<int>(element_region_ != nullptr);

    logging::error(active_regions == 1,
        "TEMPERATURE: exactly one node, surface or element region must be configured");
    logging::error(std::isfinite(temperature_),
        "TEMPERATURE: prescribed temperature must be finite");
    logging::error(model_data.positions != nullptr,
        "TEMPERATURE: model positions are not initialized");

    // Track the global node domain explicitly. `selected` provides O(1)
    // duplicate suppression while `node_ids` preserves deterministic traversal
    // order for the generated constraint rows.
    const Index       node_count = model_data.positions->rows;
    std::vector<bool> selected(node_count, false);
    std::vector<ID>   node_ids{};

    // Reuse one local insertion operation for all three target representations
    auto add_node = [&](ID node_id) {
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < node_count,
            "TEMPERATURE: node ", node_id, " is outside the compiled node domain");

        const Index node = static_cast<Index>(node_id);
        if (selected[node]) {
            return;
        }

        selected[node] = true;
        node_ids.push_back(node_id);
    };

    // A node-region target already contains the final scalar DOF locations
    if (node_region_) {
        for (ID node_id : *node_region_) {
            add_node(node_id);
        }
    }

    // Expand selected elements to their global nodal connectivity. Shared nodes
    // between adjacent elements are retained only on their first encounter.
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

    // Expand selected surfaces in the same manner. This is especially important
    // for connected surface patches, where neighboring faces commonly share
    // edge nodes.
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

    // Append one unit constraint row for every resolved thermal DOF:
    //
    //     [1] T_i = T_bar
    equations.reserve(equations.size() + node_ids.size());
    for (ID node_id : node_ids) {
        const constraint::EquationEntry entry{node_id, Dim(0), Precision(1)};
        equations.emplace_back(
            std::initializer_list<constraint::EquationEntry>{entry},
            temperature_
        );
    }
}

/**
 * Builds a compact diagnostic representation of the prescribed temperature.
 *
 * The output reports the semantic region without expanding its connectivity and
 * includes the absolute temperature stored by the condition.
 *
 * @return Human-readable temperature-condition description.
 */
std::string Temperature::str() const {
    // Resolve the configured semantic target to a compact type-qualified label
    const auto target = [&]() -> std::string {
        if (node_region_) {
            return "NSET " + node_region_->name + " (" + std::to_string(node_region_->size()) + ")";
        }
        if (surface_region_) {
            return "SFSET " + surface_region_->name + " (" + std::to_string(surface_region_->size()) + ")";
        }
        if (element_region_) {
            return "ELSET " + element_region_->name + " (" + std::to_string(element_region_->size()) + ")";
        }
        return std::string("(unknown)");
    }();

    std::ostringstream os;
    os << "TEMPERATURE: target=" << target << ", value=" << temperature_;
    return os.str();
}

} // namespace fem::bc
