/**
 * @file model.cpp
 * @brief Implements model-level registration, set population, lifecycle and compact reporting.
 *
 * This file contains the small behavioral operations that act directly on
 * `ModelData`: named node, element and surface set population from string
 * references, nonlinear step-cache coordination, registration of loads,
 * supports and amplitudes, and the compact stream representation of a model.
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
 * @date 25.08.2026
 */

#include "model.h"

#include "../bc/load_collector.h"
#include "../bc/support_collector.h"

#include <charconv>
#include <iterator>
#include <string>
#include <system_error>

namespace fem::model {

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
 * Adds one node reference to a named node set.
 *
 * The destination and source are expressed entirely through their deck-facing
 * string names. `source` may name another node set or one node identifier.
 * Before compilation the operation addresses the currently active Part and
 * stores sparse local identifiers directly. After compilation it addresses the
 * assembly node-set registry and resolves optional `instance` qualification
 * through the corresponding local-to-global node map.
 *
 * A source that already contains an `INSTANCE.` prefix keeps that explicit
 * qualification instead of receiving the separate `instance` argument. Copying
 * a set onto itself is treated as a no-op.
 *
 * @param set Destination node-set name.
 * @param source Source node-set name or node identifier.
 * @param instance Optional Instance name used for an unqualified assembly source.
 */
void Model::add_nodes_to_set(const std::string& set,
                             const std::string& source,
                             const std::string& instance) {
    // Validate the model and the destination name before selecting identifier space
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(!set.empty(),
        "Model: node set name cannot be empty");
    logging::error(_data->compiled || instance.empty(),
        "Model: node set Instance qualification requires a compiled model");

    if (!_data->compiled) {
        // Resolve the destination in the active semantic Part
        const auto part = _data->parts.get();
        logging::error(part != nullptr,
            "Model: no active part is available for node set ", set);
        auto destination = part->node_sets.activate(set);

        // Copy a named local set when the source resolves in the same Part
        if (part->node_sets.has(source)) {
            const auto source_set = part->node_sets.get(source);
            logging::error(source_set != nullptr,
                "Model: node set ", source, " is not initialized");
            if (source_set != destination) destination->add(*source_set);
            return;
        }

        // Otherwise interpret the source as one sparse local node identifier
        ID id{};
        const char* begin = source.data();
        const char* end   = begin + source.size();
        const auto [ptr, ec] = std::from_chars(begin, end, id);
        logging::error(ec == std::errc{} && ptr == end,
            "Model: node reference '", source, "' is neither a node set nor a valid node identifier");
        destination->add(id);
        return;
    }

    // Resolve optional assembly Instance qualification without overriding an
    // already qualified source such as INSTANCE.NSET or INSTANCE.ID
    logging::error(instance.empty() || _data->instances.has(instance),
        "Model: instance ", instance, " is not defined");
    const std::string qualified = instance.empty() || source.find('.') != std::string::npos
        ? source : instance + "." + source;
    auto destination = _data->node_sets.activate(set);

    // Copy an existing assembly node set when the qualified name resolves directly
    if (_data->node_sets.has(qualified)) {
        const auto source_set = _data->node_sets.get(qualified);
        logging::error(source_set != nullptr,
            "Model: node set ", qualified, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Split an entity reference into Instance and local identifier components
    std::string resolved_instance = DEFAULT_INSTANCE_NAME;
    std::string local             = qualified;
    const auto separator = qualified.rfind('.');
    if (separator != std::string::npos) {
        resolved_instance = qualified.substr(0, separator);
        local             = qualified.substr(separator + 1);
        logging::error(!resolved_instance.empty() && !local.empty(),
            "Model: node reference '", qualified, "' must use INSTANCE.ID");
    }

    // Convert the local identifier and map it into dense assembly node space
    ID id{};
    const char* begin = local.data();
    const char* end   = begin + local.size();
    const auto [ptr, ec] = std::from_chars(begin, end, id);
    logging::error(ec == std::errc{} && ptr == end,
        "Model: node reference '", qualified, "' is neither a node set nor a valid node identifier");
    destination->add(compiled_node_id(id, resolved_instance));
}

/**
 * Adds one element reference to a named element set.
 *
 * `source` may refer to another element set or to one element identifier. The
 * active Part owns the operation before compilation. Afterwards the compiled
 * assembly element-set registry is used and optional Instance qualification is
 * resolved through the corresponding element identifier map.
 *
 * A source that is already qualified as `INSTANCE.SET` or `INSTANCE.ID` keeps
 * that qualification. Copying a set onto itself leaves the destination unchanged.
 *
 * @param set Destination element-set name.
 * @param source Source element-set name or element identifier.
 * @param instance Optional Instance name used for an unqualified assembly source.
 */
void Model::add_elements_to_set(const std::string& set,
                                const std::string& source,
                                const std::string& instance) {
    // Validate the model and the destination name before selecting identifier space
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(!set.empty(),
        "Model: element set name cannot be empty");
    logging::error(_data->compiled || instance.empty(),
        "Model: element set Instance qualification requires a compiled model");

    if (!_data->compiled) {
        // Resolve the destination in the active semantic Part
        const auto part = _data->parts.get();
        logging::error(part != nullptr,
            "Model: no active part is available for element set ", set);
        auto destination = part->elem_sets.activate(set);

        // Copy a named local set when the source resolves in the same Part
        if (part->elem_sets.has(source)) {
            const auto source_set = part->elem_sets.get(source);
            logging::error(source_set != nullptr,
                "Model: element set ", source, " is not initialized");
            if (source_set != destination) destination->add(*source_set);
            return;
        }

        // Otherwise interpret the source as one sparse local element identifier
        ID id{};
        const char* begin = source.data();
        const char* end   = begin + source.size();
        const auto [ptr, ec] = std::from_chars(begin, end, id);
        logging::error(ec == std::errc{} && ptr == end,
            "Model: element reference '", source, "' is neither an element set nor a valid element identifier");
        destination->add(id);
        return;
    }

    // Resolve optional assembly Instance qualification without overriding an
    // already qualified source such as INSTANCE.ELSET or INSTANCE.ID
    logging::error(instance.empty() || _data->instances.has(instance),
        "Model: instance ", instance, " is not defined");
    const std::string qualified = instance.empty() || source.find('.') != std::string::npos
        ? source : instance + "." + source;
    auto destination = _data->elem_sets.activate(set);

    // Copy an existing assembly element set when the qualified name resolves directly
    if (_data->elem_sets.has(qualified)) {
        const auto source_set = _data->elem_sets.get(qualified);
        logging::error(source_set != nullptr,
            "Model: element set ", qualified, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Split an entity reference into Instance and local identifier components
    std::string resolved_instance = DEFAULT_INSTANCE_NAME;
    std::string local             = qualified;
    const auto separator = qualified.rfind('.');
    if (separator != std::string::npos) {
        resolved_instance = qualified.substr(0, separator);
        local             = qualified.substr(separator + 1);
        logging::error(!resolved_instance.empty() && !local.empty(),
            "Model: element reference '", qualified, "' must use INSTANCE.ID");
    }

    // Convert the local identifier and map it into dense assembly element space
    ID id{};
    const char* begin = local.data();
    const char* end   = begin + local.size();
    const auto [ptr, ec] = std::from_chars(begin, end, id);
    logging::error(ec == std::errc{} && ptr == end,
        "Model: element reference '", qualified, "' is neither an element set nor a valid element identifier");
    destination->add(compiled_element_id(id, resolved_instance));
}

/**
 * Adds one surface reference to a named surface set.
 *
 * `source` may name another surface set or one surface identifier. Before
 * compilation references remain in the active Part's sparse surface space.
 * After compilation they are added to the assembly surface-set registry and
 * optional Instance qualification is translated through the compiled surface map.
 *
 * Explicitly qualified sources keep their own Instance prefix. Copying a set
 * onto itself is ignored rather than duplicating its entries.
 *
 * @param set Destination surface-set name.
 * @param source Source surface-set name or surface identifier.
 * @param instance Optional Instance name used for an unqualified assembly source.
 */
void Model::add_surfaces_to_set(const std::string& set,
                                const std::string& source,
                                const std::string& instance) {
    // Validate the model and the destination name before selecting identifier space
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(!set.empty(),
        "Model: surface set name cannot be empty");
    logging::error(_data->compiled || instance.empty(),
        "Model: surface set Instance qualification requires a compiled model");

    if (!_data->compiled) {
        // Resolve the destination in the active semantic Part
        const auto part = _data->parts.get();
        logging::error(part != nullptr,
            "Model: no active part is available for surface set ", set);
        auto destination = part->surface_sets.activate(set);

        // Copy a named local set when the source resolves in the same Part
        if (part->surface_sets.has(source)) {
            const auto source_set = part->surface_sets.get(source);
            logging::error(source_set != nullptr,
                "Model: surface set ", source, " is not initialized");
            if (source_set != destination) destination->add(*source_set);
            return;
        }

        // Otherwise interpret the source as one sparse local surface identifier
        ID id{};
        const char* begin = source.data();
        const char* end   = begin + source.size();
        const auto [ptr, ec] = std::from_chars(begin, end, id);
        logging::error(ec == std::errc{} && ptr == end,
            "Model: surface reference '", source, "' is neither a surface set nor a valid surface identifier");
        destination->add(id);
        return;
    }

    // Resolve optional assembly Instance qualification without overriding an
    // already qualified source such as INSTANCE.SFSET or INSTANCE.ID
    logging::error(instance.empty() || _data->instances.has(instance),
        "Model: instance ", instance, " is not defined");
    const std::string qualified = instance.empty() || source.find('.') != std::string::npos
        ? source : instance + "." + source;
    auto destination = _data->surface_sets.activate(set);

    // Copy an existing assembly surface set when the qualified name resolves directly
    if (_data->surface_sets.has(qualified)) {
        const auto source_set = _data->surface_sets.get(qualified);
        logging::error(source_set != nullptr,
            "Model: surface set ", qualified, " is not initialized");
        if (source_set != destination) destination->add(*source_set);
        return;
    }

    // Split an entity reference into Instance and local identifier components
    std::string resolved_instance = DEFAULT_INSTANCE_NAME;
    std::string local             = qualified;
    const auto separator = qualified.rfind('.');
    if (separator != std::string::npos) {
        resolved_instance = qualified.substr(0, separator);
        local             = qualified.substr(separator + 1);
        logging::error(!resolved_instance.empty() && !local.empty(),
            "Model: surface reference '", qualified, "' must use INSTANCE.ID");
    }

    // Convert the local identifier and map it into dense assembly surface space
    ID id{};
    const char* begin = local.data();
    const char* end   = begin + local.size();
    const auto [ptr, ec] = std::from_chars(begin, end, id);
    logging::error(ec == std::errc{} && ptr == end,
        "Model: surface reference '", qualified, "' is neither a surface set nor a valid surface identifier");
    destination->add(compiled_surface_id(id, resolved_instance));
}

/**
 * Transfers one load into the currently active load collector.
 *
 * Loads reference compiled assembly regions, so registration is permitted only
 * after the semantic topology has been flattened. Ownership of the polymorphic
 * load is moved into the active collector.
 *
 * @param load Load definition to register.
 */
void Model::add_load(bc::Load::Ptr load) {
    // Validate the compiled model state, load ownership and collector target
    logging::error(_data != nullptr && _data->compiled,
        "Model: loads require a compiled model");
    logging::error(load != nullptr,
        "Model: cannot add a null load");
    logging::error(_data->load_cols.has_any() && _data->load_cols.get() != nullptr,
        "Model: no load collector is active");

    // Transfer the validated load into the selected collector
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
 * Transfers one support into the currently active support collector.
 *
 * Supports resolve compiled assembly regions when constraint equations are
 * collected. Registration therefore requires a compiled model and an active
 * collector. The support value is moved into that collector.
 *
 * @param support Support definition to register.
 */
void Model::add_support(bc::Support support) {
    // Validate the compiled model state and collector target
    logging::error(_data != nullptr && _data->compiled,
        "Model: supports require a compiled model");
    logging::error(_data->supp_cols.has_any() && _data->supp_cols.get() != nullptr,
        "Model: no support collector is active");

    // Transfer the support into the selected collector
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
