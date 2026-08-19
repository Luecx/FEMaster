/**
 * @file model_overview.cpp
 * @brief Implements the structured diagnostic overview of a FEM model.
 *
 * The overview exposes both sides of the model architecture: reusable semantic
 * part/instance topology and the dense assembly created by `Model::compile()`.
 * It also reports named regions, shared definitions, assignments, constraints
 * and boundary-condition collectors through the hierarchical FEMaster logger.
 *
 * The compact stream representation remains implemented separately in
 * `model.cpp` and has no logging side effects.
 *
 * @see Model
 * @see ModelData
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "model.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "../core/logging.h"

namespace fem::model {

/**
 * Logs a deterministic hierarchical overview of the model representation.
 *
 * The report first describes reusable parts and their rigid instances before
 * presenting the dense compiled assembly. Compiled elements are grouped by
 * runtime type, and every named region and shared definition is listed in
 * alphabetical order. Sections, features, constraints and load/support
 * collectors complete the report.
 *
 * Logger indentation represents ownership and grouping relationships. Null
 * entries are identified explicitly so that partially populated model state
 * remains diagnosable instead of being silently omitted.
 */
void Model::print_overview() const {
    const auto& model_data = *_data;

    // Produce stable diagnostic output even though the owning dictionaries use
    // unordered storage internally.
    auto sorted_names = [](const auto& entries) {
        std::vector<std::string> names;
        for (const auto& [name, entry] : entries) {
            (void) entry;
            names.push_back(name);
        }
        std::sort(names.begin(), names.end());
        return names;
    };

    // Count registry entries without depending on their concrete storage type
    auto entry_count = [](const auto& entries) {
        return static_cast<Index>(std::distance(entries.begin(), entries.end()));
    };

    // Print one family of compiled regions. Region sizes count dense assembly
    // identifiers, including the aggregate *ALL regions owned by each registry.
    auto print_regions = [&](const char* title, const auto& sets) {
        const auto names = sorted_names(sets);

        logging::info(true, title, " (", names.size(), ")");
        logging::up();
        for (const auto& name : names) {
            const auto& region = sets._data.at(name);
            if (region) {
                logging::info(true, name, " (", region->size(), ")");
            } else {
                logging::info(true, name, " (null)");
            }
        }
        logging::down();
    };

    // Print named definitions without depending on their concrete type. The
    // dictionary key is authoritative and remains available for null entries.
    auto print_definitions = [&](const char* title, const auto& definitions) {
        const auto names = sorted_names(definitions);

        logging::info(true, title, " (", names.size(), ")");
        logging::up();
        for (const auto& name : names) {
            logging::info(true, name, definitions.get(name) ? "" : " (null)");
        }
        logging::down();
    };

    // Start one indented report rooted at the model itself
    logging::info(true, "");
    logging::info(true, "Model overview");
    logging::up();
    logging::info(true, "Compiled: ", model_data.compiled ? "true" : "false");

    // Report reusable part definitions in their local identifier spaces
    logging::info(true, "");
    logging::info(true, "Semantic topology");
    logging::up();

    const auto part_names = sorted_names(model_data.parts);
    logging::info(true, "Parts (", part_names.size(), ")");
    logging::up();
    for (const auto& name : part_names) {
        const auto part = model_data.parts.get(name);
        logging::info(true, name, part ? "" : " (null)");
        if (!part) {
            continue;
        }

        logging::up();
        logging::info(true, "Nodes       : ", part->nodes.size());
        logging::info(true, "Elements    : ", part->elements.size());
        logging::info(true, "Surfaces    : ", part->surfaces.size());
        logging::info(true, "Lines       : ", part->lines.size());
        logging::info(true, "Sections    : ", part->sections.size());
        logging::info(true, "Node sets   : ", entry_count(part->node_sets));
        logging::info(true, "Element sets: ", entry_count(part->elem_sets));
        logging::info(true, "Surface sets: ", entry_count(part->surface_sets));
        logging::info(true, "Line sets   : ", entry_count(part->line_sets));
        logging::down();
    }
    logging::down();

    // Report instance ownership and the number of local identifiers mapped into
    // dense assembly storage by the compile pass.
    const auto instance_names = sorted_names(model_data.instances);
    logging::info(true, "Instances (", instance_names.size(), ")");
    logging::up();
    for (const auto& name : instance_names) {
        const auto instance = model_data.instances.get(name);
        logging::info(true, name, instance ? "" : " (null)");
        if (!instance) {
            continue;
        }

        logging::up();
        logging::info(true, "Part    : ", instance->part ? instance->part->name : "(null)");
        logging::info(true, "Nodes   : ", instance->node_ids.size());
        logging::info(true, "Elements: ", instance->element_ids.size());
        logging::info(true, "Surfaces: ", instance->surface_ids.size());
        logging::info(true, "Lines   : ", instance->line_ids.size());
        logging::down();
    }
    logging::down();
    logging::down();

    // Count dense elements by their runtime type. Empty solver slots remain
    // visible as UNDEFINED instead of disappearing from the diagnostic total.
    std::map<std::string, Index> element_types;
    for (const auto& element : model_data.elements) {
        std::string type = element ? element->type_name() : "";
        if (type.empty()) {
            type = "UNDEFINED";
        }
        ++element_types[type];
    }

    // Report the dense solver-facing assembly created by Model::compile()
    const Index node_count = model_data.positions ? model_data.positions->rows : Index(0);

    logging::info(true, "");
    logging::info(true, "Compiled assembly");
    logging::up();
    logging::info(true, "Nodes (", node_count, ")");
    logging::info(true, "Elements (", model_data.elements.size(), ")");
    logging::up();
    for (const auto& [type, count] : element_types) {
        logging::info(true, type, " (", count, ")");
    }
    logging::down();
    logging::info(true, "Surfaces (", model_data.surfaces.size(), ")");
    logging::info(true, "Lines (", model_data.lines.size(), ")");
    logging::down();

    // List all compiled region namespaces and their dense entity counts
    logging::info(true, "");
    logging::info(true, "Regions");
    logging::up();
    print_regions("Node sets",    model_data.node_sets);
    print_regions("Element sets", model_data.elem_sets);
    print_regions("Surface sets", model_data.surface_sets);
    print_regions("Line sets",    model_data.line_sets);
    logging::down();

    // List shared named resources independently of topology and assignment data
    logging::info(true, "");
    logging::info(true, "Definitions");
    logging::up();
    print_definitions("Materials",          model_data.materials);
    print_definitions("Profiles",           model_data.profiles);
    print_definitions("Coordinate systems", model_data.coordinate_systems);
    print_definitions("Amplitudes",         model_data.amplitudes);
    logging::down();

    // Report compiled section assignments and non-element operator features
    logging::info(true, "");
    logging::info(true, "Sections and features");
    logging::up();
    logging::info(true, "Sections (", model_data.sections.size(), ")");
    logging::up();
    for (const auto& section : model_data.sections) {
        logging::info(true, section ? section->str() : "Section: (null)");
    }
    logging::down();
    logging::info(true, "Features (", model_data.features.size(), ")");
    logging::up();
    const Index undefined_features = static_cast<Index>(std::count(
        model_data.features.begin(),
        model_data.features.end(),
        nullptr
    ));
    logging::info(undefined_features > 0, "Undefined: ", undefined_features);
    logging::down();
    logging::down();

    // Summarize every stored constraint family. Couplings expose a stable
    // printable representation; the remaining types currently expose counts.
    logging::info(true, "");
    logging::info(true, "Constraints");
    logging::up();
    logging::info(true, "Connectors: ", model_data.connectors.size());
    logging::info(true, "Couplings : ", model_data.couplings.size());
    logging::up();
    for (const auto& coupling : model_data.couplings) {
        logging::info(true, coupling.str());
    }
    logging::down();
    logging::info(true, "Ties      : ", model_data.ties.size());
    logging::info(true, "Contacts  : ", model_data.contacts.size());
    logging::info(true, "RBMs      : ", model_data.rbms.size());
    logging::info(true, "Equations : ", model_data.equations.size());
    logging::down();

    // Expand support collectors and their value entries in deterministic name
    // order. Supports are stored by value and are therefore never null.
    const auto support_collector_names = sorted_names(model_data.supp_cols);

    logging::info(true, "");
    logging::info(true, "Support collectors (", support_collector_names.size(), ")");
    logging::up();
    for (const auto& name : support_collector_names) {
        const auto& collector = model_data.supp_cols._data.at(name);
        if (!collector) {
            logging::info(true, name, " (null)");
            continue;
        }

        logging::info(true, name, " (", collector->size(), ")");
        logging::up();
        for (const auto& support : collector->entries()) {
            logging::info(true, support.str());
        }
        logging::down();
    }
    logging::down();

    // Expand load collectors and identify missing polymorphic load entries
    // explicitly instead of silently skipping them.
    const auto load_collector_names = sorted_names(model_data.load_cols);

    logging::info(true, "");
    logging::info(true, "Load collectors (", load_collector_names.size(), ")");
    logging::up();
    for (const auto& name : load_collector_names) {
        const auto& collector = model_data.load_cols._data.at(name);
        if (!collector) {
            logging::info(true, name, " (null)");
            continue;
        }

        logging::info(true, name, " (", collector->size(), ")");
        logging::up();
        for (const auto& load : collector->entries()) {
            logging::info(true, load ? load->str() : "Load: (null)");
        }
        logging::down();
    }
    logging::down();

    // Close the model-level indentation established at the start of the report
    logging::down();
}

} // namespace fem::model
