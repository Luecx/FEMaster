/**
 * @file model.cpp
 * @brief Implements model-level registration, set population, lifecycle and compact reporting.
 *
 * This file contains the small behavioral operations that act directly on
 * `ModelData`: insertion of named-set or scalar references into existing
 * regions, nonlinear step-cache coordination, registration of loads, supports
 * and amplitudes, and the compact stream representation of a model.
 *
 * Named set population remains independent of the parser format. Before
 * compilation references address the active Part and retain sparse local IDs;
 * afterwards they may use optional Instance qualification and are resolved
 * through the compiled local-to-global maps.
 *
 * Part/instance flattening remains in `model_compile.cpp`, global matrix and
 * load construction remain in `model_build.cpp`, and the detailed hierarchical
 * diagnostic report is implemented in `model_overview.cpp`.
 *
 * @see Model
 * @see ModelData
 * @see Part
 * @see Instance
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "model.h"

#include "../bc/load_collector.h"
#include "../bc/support_collector.h"

#include <charconv>
#include <iterator>
#include <string>
#include <system_error>

namespace fem::model {

namespace {

/**
 * Parses one sparse Part-local entity identifier.
 *
 * The complete source token must represent an integer. Named-region lookup is
 * performed by the calling set operation before this scalar parser is used.
 *
 * @param source Local identifier token.
 * @param entity Entity name used in validation diagnostics.
 * @return Parsed sparse local identifier.
 */
ID parse_local_id(const std::string& source, const char* entity) {
    // Convert the complete token without accepting trailing characters
    ID id{};
    const char* begin = source.data();
    const char* end   = begin + source.size();
    const auto [ptr, ec] = std::from_chars(begin, end, id);

    logging::error(ec == std::errc{} && ptr == end,
        "Model: ", entity, " reference '", source, "' is not a valid local identifier");
    return id;
}

} // namespace

/**
 * Initializes persistent element state required during one analysis step.
 *
 * Every compiled structural element receives a lifecycle callback. Elements
 * that do not require persistent step data use the default no-op implementation,
 * while formulations such as finite-rotation shells construct their cached
 * reference geometry here.
 *
 * Initialization is transactional at model level: if one element throws, the
 * method invokes `step_end()` for the complete model to release every cache that
 * may already have been created, then propagates the original exception.
 */
void Model::step_begin() {
    // Initialize analysis-local state for every compiled structural element
    try {
        for (auto& elem : _data->elements) {
            if (elem) {
                if (auto* structural = elem->as<StructuralElement>()) {
                    structural->step_begin();
                }
            }
        }
    } catch (...) {
        // Release partially initialized state before propagating the failure
        step_end();
        throw;
    }
}

/**
 * Releases persistent element state associated with the current analysis step.
 *
 * The callback is sent to every compiled structural element and is safe after a
 * partially completed `step_begin()`. Persistent model topology, material
 * history and reusable thread-local evaluation workspaces remain unchanged.
 */
void Model::step_end() {
    // Release only the analysis-local caches owned by structural elements
    for (auto& elem : _data->elements) {
        if (elem) {
            if (auto* structural = elem->as<StructuralElement>()) {
                structural->step_end();
            }
        }
    }
}

/**
 * Adds one local node-set or node identifier to a set in the active Part.
 *
 * A named source set is copied within the same Part namespace. Otherwise the
 * source is parsed as one sparse local node identifier and inserted directly.
 *
 * @param set Destination node-set name in the active Part.
 * @param source Local node-set name or node identifier.
 */
void Model::add_nodes_to_part_set(const std::string& set, const std::string& source) {
    // Validate the semantic identifier space and active Part
    logging::error(_data != nullptr && !_data->compiled,
        "Model: part node sets cannot be modified after compile()");
    logging::error(!set.empty(),
        "Model: node set name cannot be empty");

    const auto part = _data->parts.get();
    logging::error(part != nullptr,
        "Model: no active part is available for node set ", set);

    // Select the destination and prefer a named source set over scalar parsing
    const auto destination = part->node_sets.activate(set);
    if (part->node_sets.has(source)) {
        const auto source_set = part->node_sets.get(source);
        logging::error(source_set != nullptr,
            "Model: node set ", source, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Retain one sparse Part-local node identifier
    destination->add(parse_local_id(source, "node"));
}

/**
 * Adds one compiled node-set or scalar node reference to an assembly set.
 *
 * A named source set is copied within assembly node space. Otherwise `source`
 * is resolved as `ID` or `INSTANCE.ID` through the compiled node maps.
 *
 * @param set Destination assembly node-set name.
 * @param source Compiled node-set name or scalar node reference.
 */
void Model::add_nodes_to_assembly_set(const std::string& set, const std::string& source) {
    // Validate the compiled identifier space and destination name
    logging::error(_data != nullptr && _data->compiled,
        "Model: assembly node sets require a compiled model");
    logging::error(!set.empty(),
        "Model: node set name cannot be empty");

    // Select the destination and prefer a named compiled source set
    const auto destination = _data->node_sets.activate(set);
    if (_data->node_sets.has(source)) {
        const auto source_set = _data->node_sets.get(source);
        logging::error(source_set != nullptr,
            "Model: node set ", source, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Map one semantic node reference into dense assembly space
    destination->add(compiled_node_id(source));
}

/**
 * Adds one local element-set or element identifier to a set in the active Part.
 *
 * A named source set is copied within the same Part namespace. Otherwise the
 * source is parsed as one sparse local element identifier and inserted directly.
 *
 * @param set Destination element-set name in the active Part.
 * @param source Local element-set name or element identifier.
 */
void Model::add_elements_to_part_set(const std::string& set, const std::string& source) {
    // Validate the semantic identifier space and active Part
    logging::error(_data != nullptr && !_data->compiled,
        "Model: part element sets cannot be modified after compile()");
    logging::error(!set.empty(),
        "Model: element set name cannot be empty");

    const auto part = _data->parts.get();
    logging::error(part != nullptr,
        "Model: no active part is available for element set ", set);

    // Select the destination and prefer a named source set over scalar parsing
    const auto destination = part->elem_sets.activate(set);
    if (part->elem_sets.has(source)) {
        const auto source_set = part->elem_sets.get(source);
        logging::error(source_set != nullptr,
            "Model: element set ", source, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Retain one sparse Part-local element identifier
    destination->add(parse_local_id(source, "element"));
}

/**
 * Adds one compiled element-set or scalar element reference to an assembly set.
 *
 * A named source set is copied within assembly element space. Otherwise
 * `source` is resolved as `ID` or `INSTANCE.ID` through the compiled element
 * maps.
 *
 * @param set Destination assembly element-set name.
 * @param source Compiled element-set name or scalar element reference.
 */
void Model::add_elements_to_assembly_set(const std::string& set, const std::string& source) {
    // Validate the compiled identifier space and destination name
    logging::error(_data != nullptr && _data->compiled,
        "Model: assembly element sets require a compiled model");
    logging::error(!set.empty(),
        "Model: element set name cannot be empty");

    // Select the destination and prefer a named compiled source set
    const auto destination = _data->elem_sets.activate(set);
    if (_data->elem_sets.has(source)) {
        const auto source_set = _data->elem_sets.get(source);
        logging::error(source_set != nullptr,
            "Model: element set ", source, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Map one semantic element reference into dense assembly space
    destination->add(compiled_element_id(source));
}

/**
 * Adds one local surface-set or surface identifier to a set in the active Part.
 *
 * A named source set is copied within the same Part namespace. Otherwise the
 * source is parsed as one sparse local surface identifier and inserted directly.
 *
 * @param set Destination surface-set name in the active Part.
 * @param source Local surface-set name or surface identifier.
 */
void Model::add_surfaces_to_part_set(const std::string& set, const std::string& source) {
    // Validate the semantic identifier space and active Part
    logging::error(_data != nullptr && !_data->compiled,
        "Model: part surface sets cannot be modified after compile()");
    logging::error(!set.empty(),
        "Model: surface set name cannot be empty");

    const auto part = _data->parts.get();
    logging::error(part != nullptr,
        "Model: no active part is available for surface set ", set);

    // Select the destination and prefer a named source set over scalar parsing
    const auto destination = part->surface_sets.activate(set);
    if (part->surface_sets.has(source)) {
        const auto source_set = part->surface_sets.get(source);
        logging::error(source_set != nullptr,
            "Model: surface set ", source, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Retain one sparse Part-local surface identifier
    destination->add(parse_local_id(source, "surface"));
}

/**
 * Adds one compiled surface-set or scalar surface reference to an assembly set.
 *
 * A named source set is copied within assembly surface space. Otherwise
 * `source` is resolved as `ID` or `INSTANCE.ID` through the compiled surface
 * maps.
 *
 * @param set Destination assembly surface-set name.
 * @param source Compiled surface-set name or scalar surface reference.
 */
void Model::add_surfaces_to_assembly_set(const std::string& set, const std::string& source) {
    // Validate the compiled identifier space and destination name
    logging::error(_data != nullptr && _data->compiled,
        "Model: assembly surface sets require a compiled model");
    logging::error(!set.empty(),
        "Model: surface set name cannot be empty");

    // Select the destination and prefer a named compiled source set
    const auto destination = _data->surface_sets.activate(set);
    if (_data->surface_sets.has(source)) {
        const auto source_set = _data->surface_sets.get(source);
        logging::error(source_set != nullptr,
            "Model: surface set ", source, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Map one semantic surface reference into dense assembly space
    destination->add(compiled_surface_id(source));
}

/**
 * Transfers one Neumann condition into the currently active load collector.
 *
 * Natural boundary conditions reference compiled assembly regions, so
 * registration is permitted only after the semantic topology has been flattened.
 * Ownership of the polymorphic condition is moved into the active collector.
 *
 * @param load Neumann condition to register.
 */
void Model::add_load(bc::Neumann::Ptr load) {
    // Validate the compiled model state, condition ownership and collector target
    logging::error(_data != nullptr && _data->compiled,
        "Model: loads require a compiled model");
    logging::error(load != nullptr,
        "Model: cannot add a null load");
    logging::error(_data->load_cols.has_any() && _data->load_cols.get() != nullptr,
        "Model: no load collector is active");

    // Transfer the validated condition into the selected collector
    _data->load_cols.get()->add(std::move(load));
}

/**
 * Registers one named amplitude as a shared model definition.
 *
 * Amplitudes are independent of part/instance identifier flattening and may be
 * added before or after compilation. Names are unique within the model-wide
 * amplitude dictionary, which takes shared ownership of the supplied object.
 *
 * @param amplitude Named time-history definition to register.
 */
void Model::add_amplitude(bc::Amplitude::Ptr amplitude) {
    // Validate ownership and the model-wide amplitude namespace
    logging::error(amplitude != nullptr,
        "Model: cannot add a null amplitude");
    logging::error(!_data->amplitudes.has(amplitude->name),
        "Model: amplitude ", amplitude->name, " is already defined");

    // Transfer the unique named definition into persistent model storage
    _data->amplitudes.add(std::move(amplitude));
}

/**
 * Transfers one Dirichlet condition into the currently active support collector.
 *
 * Essential conditions resolve compiled assembly regions when equations are
 * collected. Registration therefore requires a compiled model and an active
 * collector. Structural supports and prescribed temperatures share this
 * polymorphic ownership path.
 *
 * @param support Dirichlet condition to register.
 */
void Model::add_support(bc::Dirichlet::Ptr support) {
    // Validate the compiled model state, condition ownership and collector target
    logging::error(_data != nullptr && _data->compiled,
        "Model: supports require a compiled model");
    logging::error(support != nullptr,
        "Model: cannot add a null support");
    logging::error(_data->supp_cols.has_any() && _data->supp_cols.get() != nullptr,
        "Model: no support collector is active");

    // Transfer the condition into the selected collector
    _data->supp_cols.get()->add(std::move(support));
}

/**
 * Writes a compact model summary to the supplied stream.
 *
 * The representation contains only stable high-level counts and the topology
 * compilation state. Detailed hierarchical diagnostics remain the
 * responsibility of `Model::print_overview()` and are not emitted through the
 * logger as a side effect of streaming.
 *
 * @param ostream Destination stream receiving the model summary.
 * @param model Model whose semantic and compiled entity counts are reported.
 * @return Reference to `ostream` for chained stream operations.
 */
std::ostream& operator<<(std::ostream& ostream, const Model& model) {
    const auto& model_data = *model._data;

    // Count named semantic containers through their public iterator interface.
    // This avoids depending on the associative storage selected by Dict.
    auto entry_count = [](const auto& entries) {
        return static_cast<Index>(std::distance(entries.begin(), entries.end()));
    };

    // Use the compiled position field as the authoritative dense node count
    const Index node_count = model_data.positions ? model_data.positions->rows : Index(0);

    // Write a compact and side-effect-free summary with aligned labels
    ostream << "compiled  = " << (model_data.compiled ? "true" : "false") << '\n';
    ostream << "parts     = " << entry_count(model_data.parts)            << '\n';
    ostream << "instances = " << entry_count(model_data.instances)        << '\n';
    ostream << "nodes     = " << node_count                               << '\n';
    ostream << "elements  = " << model_data.elements.size()               << '\n';
    ostream << "surfaces  = " << model_data.surfaces.size()               << '\n';
    ostream << "lines     = " << model_data.lines.size()                  << '\n';
    return ostream;
}

} // namespace fem::model
