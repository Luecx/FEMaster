/**
 * @file model_compile.cpp
 * @brief Compiles part/instance topology and local assignments into solver storage.
 *
 * `Model::compile()` is intentionally narrow: it freezes only the semantic
 * part/instance data that must be flattened into the assembly. Nodes, elements,
 * surfaces and lines are copied from every instance, transformed and rewired to
 * dense global identifiers. Instance-local sets and section regions are
 * materialized in the same global identifier space.
 *
 * Definition objects such as materials, profiles and coordinate systems are not
 * part of that one-way topology transition. They may be registered before or
 * after compilation. Sections may likewise be added after compilation when they
 * already reference a compiled element set; pre-compile sections remain part
 * definitions and are copied once per instance here.
 *
 * Section orientations require the same rigid placement as their owning part.
 * A compiled section therefore receives an independent transformed coordinate
 * system: orientation axes rotate with the instance and spatial origins also
 * receive the instance translation. The type-specific copy is performed by
 * `Section::instantiate()`.
 *
 * Dense identifiers are deterministic: instances are enumerated by instance id
 * and every part-local entity identifier is sorted before enumeration. The
 * resulting Instance maps remain available for resolving qualified local
 * references after the model has crossed the one-way compilation boundary.
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

#include <algorithm>
#include <charconv>
#include <cmath>
#include <limits>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

namespace fem::model {

namespace {

/**
 * Parses a compact reference to one semantic entity in an Instance namespace.
 *
 * A plain local identifier is associated with the implicit identity Instance.
 * A qualified reference is split at its final separator into `INSTANCE.ID`.
 * The complete identifier token must represent one valid integer value.
 *
 * @param reference Local `ID` or qualified `INSTANCE.ID` reference.
 * @param entity Entity name used in validation diagnostics.
 * @return Instance name and local semantic identifier.
 */
std::pair<std::string, Index> split_compiled_reference(const std::string& reference, const char* entity) {
    // Separate an optional Instance prefix from the local identifier token
    std::string instance = Model::DEFAULT_INSTANCE_NAME;
    std::string local    = reference;

    const auto separator = reference.rfind('.');
    if (separator != std::string::npos) {
        instance = reference.substr(0, separator);
        local    = reference.substr(separator + 1);

        logging::error(!instance.empty() && !local.empty(),
            "Model: ", entity, " reference '", reference, "' must use INSTANCE.ID");
    }

    // Convert the complete local token without accepting trailing characters
    Index id{};
    const char* begin = local.data();
    const char* end   = begin + local.size();
    const auto [ptr, ec] = std::from_chars(begin, end, id);

    logging::error(ec == std::errc{} && ptr == end,
        "Model: ", entity, " reference '", reference, "' does not contain a valid local identifier");
    return {std::move(instance), id};
}

} // namespace

/**
 * Constructs an empty model with an implicit root part and identity instance.
 *
 * Unqualified topology uses these two semantic objects, allowing root-level and
 * explicitly instanced input to share the same compilation and identifier
 * resolution path without a separate assembly topology type.
 */
Model::Model()
    : _data(std::make_shared<ModelData>()) {
    // Unqualified deck entities live in the default part. The corresponding
    // default instance is its identity embedding into the final assembly.
    const auto default_part = _data->parts.activate(DEFAULT_PART_NAME);
    _data->instances.activate(DEFAULT_INSTANCE_NAME, default_part, Index{0});
}

/**
 * Creates and activates a named reusable part definition.
 *
 * Parts participate in identifier flattening and may therefore be added only
 * before compilation. Names must be non-empty and unique in the semantic part
 * dictionary. The new part becomes active for subsequent local topology and
 * section construction.
 *
 * @param name Unique semantic part name.
 */
void Model::add_part(const std::string& name) {
    // Parts participate directly in flattening and may therefore only be
    // created before the one allowed compile pass.
    logging::error(_data != nullptr && !_data->compiled,
        "Model: parts cannot be added after compile()");
    logging::error(!name.empty(),
        "Model: part name cannot be empty");
    logging::error(!_data->parts.has(name),
        "Model: part ", name, " is already defined");

    // Create the part and select it as the current topology destination
    _data->parts.activate(name);
}

/**
 * Adds one rigid placement of an existing semantic part.
 *
 * The instance stores its source Part and the local-to-global transformation
 * `x = R X + t`. Finiteness, orthonormality and proper-rotation requirements are
 * validated together for all instances when compilation begins. Instance names
 * are unique and remain the namespace qualifier for inherited regions and local
 * identifier resolution.
 *
 * @param name Unique instance name.
 * @param part_name Name of the existing reusable part definition.
 * @param translation Global translation applied after rotation.
 * @param rotation Proper rotation from part-local to global coordinates.
 */
void Model::add_instance(const std::string& name,
                         const std::string& part_name,
                         Vec3 translation,
                         Mat3 rotation) {
    // Instance placement determines all dense assembly identifiers and cannot
    // change after the model has been flattened.
    logging::error(_data != nullptr && !_data->compiled,
        "Model: instances cannot be added after compile()");
    logging::error(!name.empty(),
        "Model: instance name cannot be empty");
    logging::error(_data->parts.has(part_name),
        "Model: part ", part_name, " is not defined for instance ", name);
    logging::error(!_data->instances.has(name),
        "Model: instance ", name, " is already defined");

    // Store the semantic placement without duplicating its referenced topology.
    // Instance ids follow creation order after the default instance (id 0).
    _data->instances.add(std::make_shared<Instance>(
        name,
        _data->parts.get(part_name),
        static_cast<Index>(_data->instances.size()),
        translation,
        rotation
    ));
}

/**
 * Registers a unique named coordinate-system definition.
 *
 * Coordinate systems are shared resources rather than flattened topology and
 * may be added before or after compilation. The model dictionary takes shared
 * ownership of the supplied polymorphic definition.
 *
 * @param coordinate_system Coordinate system to register.
 */
void Model::add_coordinate_system(cos::CoordinateSystem::Ptr coordinate_system) {
    // Coordinate systems are definitions rather than assembly topology. Abaqus
    // *TRANSFORM deliberately creates nodal systems after compile(), so no
    // compiled-state guard belongs here.
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(coordinate_system != nullptr,
        "Model: cannot add a null coordinate system");
    logging::error(!_data->coordinate_systems.has(coordinate_system->name),
        "Model: coordinate system ", coordinate_system->name, " is already defined");

    // Transfer the unique named definition into persistent model storage
    _data->coordinate_systems.add(std::move(coordinate_system));
}

/**
 * Registers a unique named material definition.
 *
 * Materials are shared by section assignments and do not participate directly
 * in part/instance identifier flattening, so registration is independent of the
 * compile state.
 *
 * @param material Material definition to register.
 */
void Model::add_material(material::Material::Ptr material) {
    // Materials do not affect part/instance identifier flattening. They may be
    // introduced whenever a later section or analysis object needs them.
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(material != nullptr,
        "Model: cannot add a null material");
    logging::error(!_data->materials.has(material->name),
        "Model: material ", material->name, " is already defined");

    // Transfer the unique named definition into persistent model storage
    _data->materials.add(std::move(material));
}

/**
 * Registers a unique named beam-profile definition.
 *
 * Profiles are shared resources referenced by beam sections and remain
 * independent of semantic topology compilation.
 *
 * @param profile Profile definition to register.
 */
void Model::add_profile(Profile::Ptr profile) {
    // Profiles behave like materials: they are shared definitions and are not
    // part of the topology that compile() freezes.
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(profile != nullptr,
        "Model: cannot add a null profile");
    logging::error(!_data->profiles.has(profile->name),
        "Model: profile ", profile->name, " is already defined");

    // Transfer the unique named definition into persistent model storage
    _data->profiles.add(std::move(profile));
}

/**
 * Registers a section assignment in semantic or compiled element space.
 *
 * Before compilation the section region must be the exact region object owned
 * by the active Part. The section remains a semantic definition and is copied
 * once for every Instance during `compile()`. After compilation the region must
 * instead be the exact registered assembly element set; the section is stored
 * directly and immediately assigned to its target elements.
 *
 * This pointer-identity requirement prevents stale or same-named regions from a
 * different identifier space from crossing the compilation boundary.
 *
 * @param section Section assignment to register.
 */
void Model::add_section(Section::Ptr section) {
    // Validate section ownership and its required target region
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(section != nullptr,
        "Model: cannot add a null section");
    logging::error(section->region_ != nullptr,
        "Model: section has no element region");

    if (_data->compiled) {
        // Post-compile sections already live in assembly space. Requiring the
        // exact registered region pointer prevents accidentally feeding a
        // stale part-local ELSET into dense solver storage.
        logging::error(_data->elem_sets.has(section->region_->name) && _data->elem_sets.get(section->region_->name) == section->region_,
            "Model: post-compile section region ", section->region_->name,
            " is not a compiled element set");

        // Store and bind the assembly-space assignment immediately
        _data->sections.push_back(std::move(section));
        assign_sections();
        return;
    }

    // Before compilation a section belongs to the currently active part. Its
    // region must therefore be one of that part's own element sets. compile()
    // later creates one global section assignment for every instance.
    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");
    logging::error(active->elem_sets.has(section->region_->name) && active->elem_sets.get(section->region_->name) == section->region_,
        "Model: section region ", section->region_->name,
        " does not belong to active part ", active->name);

    // Retain the semantic assignment for instance-specific compilation
    active->sections.push_back(std::move(section));
}

/**
 * Flattens all semantic Parts and Instances into dense assembly storage.
 *
 * Compilation validates every rigid instance, determines exact assembly sizes
 * and allocates dense node, element, surface and line storage. Instances are
 * processed by instance id; every local entity family is enumerated by sorted
 * identifier. Nodes are transformed by `x = R X + t`, polymorphic
 * topology objects are copied and rewired through instance maps, and inherited
 * regions receive qualified assembly names.
 *
 * Part-local sections become independent assembly assignments per Instance.
 * Their element regions are remapped; spatial section data is additionally
 * transformed with the rigid instance placement where required. Part-local
 * constraints are copied the same way with remapped nodes and regions.
 *
 * Finally, sections are bound to dense elements, element-nodal/IP/MP offsets are
 * initialized and the semantic topology is frozen. Recompilation is unsupported
 * because it would invalidate every assembly-level field, region and reference.
 */
void Model::compile() {
    // Compiling twice would duplicate/renumber the assembly and invalidate every
    // global reference created after the first pass. Recompilation is therefore
    // explicitly unsupported instead of attempting cache invalidation.
    logging::error(_data != nullptr && !_data->compiled,
        "Model: compile() may only be called once");

    // Collect deterministic instance order and exact dense assembly dimensions
    std::vector<Instance::Ptr> instances;
    instances.reserve(_data->instances._data.size());

    Index total_nodes    = 0;
    Index total_elements = 0;
    Index total_surfaces = 0;
    Index total_lines    = 0;

    // Only proper rigid rotations are accepted. Reflections would change element
    // orientation and surface handedness and are not valid instance placements.
    const Precision rotation_tolerance =
        Precision(100) * std::sqrt(std::numeric_limits<Precision>::epsilon());

    // Validate every rigid placement and accumulate the dense assembly size
    for (const auto& [name, instance] : _data->instances) {
        logging::error(instance != nullptr && instance->part != nullptr,
            "Model: instance ", name, " does not reference a part");
        logging::error(instance->translation.allFinite(),
            "Model: instance ", name, " translation contains invalid values");
        logging::error(instance->rotation.allFinite(),
            "Model: instance ", name, " rotation contains invalid values");
        logging::error((instance->rotation.transpose() * instance->rotation - Mat3::Identity()).norm() <= rotation_tolerance,
            "Model: instance ", name, " rotation is not orthonormal");
        logging::error(std::abs(instance->rotation.determinant() - Precision(1)) <= rotation_tolerance,
            "Model: instance ", name, " rotation must have determinant +1");

        total_nodes    += static_cast<Index>(instance->part->nodes   .size());
        total_elements += static_cast<Index>(instance->part->elements.size());
        total_surfaces += static_cast<Index>(instance->part->surfaces.size());
        total_lines    += static_cast<Index>(instance->part->lines   .size());
        instances.push_back(instance);
    }

    // Enumerate instances by id so that the default instance (id 0) comes first
    // and the dense numbering follows the input order of the instances.
    std::sort(instances.begin(), instances.end(),
        [](const Instance::Ptr& a, const Instance::Ptr& b) {
            return a->instance_id < b->instance_id;
        });

    // Dense solver arrays are allocated exactly once from the semantic assembly
    // size. Sparse/local identifiers never escape into these arrays.
    _data->elements.assign(total_elements, nullptr);
    _data->surfaces.assign(total_surfaces, nullptr);
    _data->lines   .assign(total_lines, nullptr);
    _data->sections.clear();

    // Dense-to-semantic maps used by writers and readers to recover instance + local id
    _data->node_mapping   .resize(total_nodes);
    _data->element_mapping.resize(total_elements);

    // Rebuild the compiled set namespaces. Each Sets container automatically
    // propagates additions into its global *ALL parent collection.
    _data->node_sets    = Sets<NodeRegion>{SET_NODE_ALL};
    _data->elem_sets    = Sets<ElementRegion>{SET_ELEM_ALL};
    _data->surface_sets = Sets<SurfaceRegion>{SET_SURF_ALL};
    _data->line_sets    = Sets<LineRegion>{SET_LINE_ALL};

    // POSITION and POSITION_REFERENCE are dense nodal assembly fields. Rotations
    // occupy the remaining components later during nonlinear analyses and start
    // at zero in both fields.
    _data->fields.clear();
    _data->positions           = std::make_shared<Field>("POSITION"          , FieldDomain::NODE, total_nodes, 6);
    _data->positions_reference = std::make_shared<Field>("POSITION_REFERENCE", FieldDomain::NODE, total_nodes, 6);
    _data->positions          ->set_zero();
    _data->positions_reference->set_zero();
    _data->fields.emplace(_data->positions->name          , _data->positions);
    _data->fields.emplace(_data->positions_reference->name, _data->positions_reference);

    Index next_node    = 0;
    Index next_element = 0;
    Index next_surface = 0;
    Index next_line    = 0;

    // Part-local identifiers of one entity family in ascending order. Sorting
    // makes the dense numbering independent of the unordered_map storage in Part.
    auto sorted_ids = [](const auto& entities) {
        std::vector<Index> ids;
        ids.reserve(entities.size());
        for (const auto& entry : entities) {
            ids.push_back(entry.first);
        }
        std::sort(ids.begin(), ids.end());
        return ids;
    };

    // Materialize every instance in instance-id order
    for (const auto& instance : instances) {
        const auto& source     = instance->part;
        const auto& name       = instance->name;
        const bool  is_default = instance->instance_id == 0;

        // The default instance owns the unqualified assembly namespace. Explicit
        // instances prefix inherited part sets with "INSTANCE.".
        auto qualified_name = [&](const std::string& local_name) {
            return is_default ? local_name : name + "." + local_name;
        };

        // A section orientation belongs to the part reference configuration and
        // must follow the same rigid placement as the instance geometry. The
        // transformed coordinate system is kept section-local; the source global
        // definition remains unchanged and reusable by other instances.
        auto transformed_orientation = [&](const cos::CoordinateSystem::Ptr& orientation) {
            return orientation
                ? orientation->transformed(instance->rotation, instance->translation)
                : cos::CoordinateSystem::Ptr{};
        };

        auto& node_map    = instance->node_ids;
        auto& element_map = instance->element_ids;
        auto& surface_map = instance->surface_ids;
        auto& line_map    = instance->line_ids;

        // Copies one part-local element, surface or line family into dense
        // storage. Every copy is rewired through the node map and registered in
        // the dense *ALL region and the local-to-dense map. prepare() performs
        // the family-specific setup of the copy and returns its node count.
        auto compile_entities = [&](const auto& entities,
                                    auto&       target,
                                    const auto& all,
                                    auto&       id_map,
                                    Index&      next_id,
                                    const char* entity,
                                    auto&&      prepare) {
            for (const Index local_id : sorted_ids(entities)) {
                const Index global_id = next_id++;
                auto        copy      = entities.at(local_id)->copy();
                const Index n_nodes   = prepare(copy, global_id, local_id);

                for (Index local = 0; local < n_nodes; ++local) {
                    const Index local_node = copy->nodes()[local];
                    const auto  it         = node_map.find(local_node);
                    logging::error(it != node_map.end(),
                        "Model: ", entity, " ", local_id, " in part ", source->name,
                        " references undefined node ", local_node);
                    copy->nodes()[local] = it->second;
                }

                target[global_id] = std::move(copy);
                all->add(global_id);
                id_map.emplace(local_id, global_id);
            }
        };

        // Materializes inherited part-local regions as qualified dense regions.
        // The default *ALL regions are already global parents and must not be
        // replaced by a second same-named default-part collection.
        auto compile_sets = [&](const auto&        local_sets,
                                auto&              target_sets,
                                const std::string& all_name,
                                const auto&        id_map,
                                const char*        entity) {
            for (const auto& [local_name, local_set] : local_sets) {
                if (is_default && local_name == all_name) continue;

                auto target = target_sets.activate(qualified_name(local_name));
                for (const Index local_id : *local_set) {
                    const auto it = id_map.find(local_id);
                    logging::error(it != id_map.end(),
                        "Model: ", entity, " set ", local_name, " in part ", source->name,
                        " references undefined ", entity, " ", local_id);
                    target->add(it->second);
                }
            }
        };

        // Transform local nodal coordinates and establish the dense node map
        for (const Index local_id : sorted_ids(source->nodes)) {
            const Index global_id   = next_node++;
            const Vec3  position    = instance->rotation * source->nodes.at(local_id) + instance->translation;

            node_map.emplace(local_id, global_id);
            _data->node_mapping[global_id] = {instance, local_id};

            for (Index d = 0; d < 3; ++d) {
                (*_data->positions          )(global_id, d) = position(d);
                (*_data->positions_reference)(global_id, d) = position(d);
            }
            _data->node_sets.all()->add(global_id);
        }

        // copy() preserves the concrete element type and topology-specific
        // configuration. Assembly-owned ids, offsets, section assignment and
        // ModelData binding are always rebuilt explicitly here.
        compile_entities(source->elements, _data->elements, _data->elem_sets.all(), element_map, next_element, "element",
            [&](ElementPtr& element, Index global_id, Index local_id) {
                element->elem_id           = global_id;
                element->elem_nodal_offset = 0;
                element->elem_ip_offset    = 0;
                element->elem_mp_offset    = 0;
                element->_section          = nullptr;
                element->_model_data       = _data.get();

                _data->element_mapping[global_id] = {instance, local_id};
                return element->n_nodes();
            });

        compile_sets(source->node_sets, _data->node_sets, SET_NODE_ALL, node_map   , "node");
        compile_sets(source->elem_sets, _data->elem_sets, SET_ELEM_ALL, element_map, "element");

        compile_entities(source->surfaces, _data->surfaces, _data->surface_sets.all(), surface_map, next_surface, "surface",
            [](SurfacePtr& surface, Index, Index) { return surface->n_nodes; });
        compile_sets(source->surface_sets, _data->surface_sets, SET_SURF_ALL, surface_map, "surface");

        compile_entities(source->lines, _data->lines, _data->line_sets.all(), line_map, next_line, "line",
            [](LinePtr& line, Index, Index) { return line->n_nodes; });
        compile_sets(source->line_sets, _data->line_sets, SET_LINE_ALL, line_map, "line");

        // Sections are definitions on the source part but assignments in the
        // compiled assembly. Every instance receives its own copy with a remapped
        // element region and rigidly placed spatial section data.
        for (const Section::Ptr& source_section : source->sections) {
            auto region = std::make_shared<ElementRegion>(qualified_name(source_section->region_->name));

            for (const Index local_id : *source_section->region_) {
                const auto it = element_map.find(local_id);
                logging::error(it != element_map.end(),
                    "Model: section region ", source_section->region_->name, " in part ", source->name,
                    " references undefined element ", local_id);
                region->add(it->second);
            }

            _data->sections.push_back(source_section->instantiate(std::move(region), instance->rotation, instance->translation));
        }

        // Part-local constraints are copied once per instance. Nodes are mapped
        // through the dense node map and regions are replaced by their qualified
        // compiled counterparts materialized above.
        auto compiled_node = [&](Index local_id) {
            const auto it = node_map.find(local_id);
            logging::error(it != node_map.end(),
                "Model: constraint in part ", source->name, " references undefined node ", local_id);
            return it->second;
        };

        auto compiled_region = [&](auto& sets, const auto& region) -> std::decay_t<decltype(region)> {
            if (!region) return nullptr;
            logging::error(sets.has(qualified_name(region->name)),
                "Model: constraint region ", region->name, " in part ", source->name, " is not a named set");
            return sets.get(qualified_name(region->name));
        };

        for (auto connector : source->connectors) {
            connector.node_1_id_         = compiled_node(connector.node_1_id_);
            connector.node_2_id_         = compiled_node(connector.node_2_id_);
            connector.coordinate_system_ = transformed_orientation(connector.coordinate_system_);
            _data->connectors.push_back(std::move(connector));
        }

        for (auto mpc : source->mpcs) {
            mpc.node_1_id_ = compiled_node(mpc.node_1_id_);
            mpc.node_2_id_ = compiled_node(mpc.node_2_id_);
            _data->mpcs.push_back(std::move(mpc));
        }

        for (auto coupling : source->couplings) {
            coupling.master_node    = compiled_node(coupling.master_node);
            coupling.master_region  = compiled_region(_data->node_sets   , coupling.master_region);
            coupling.slave_nodes    = compiled_region(_data->node_sets   , coupling.slave_nodes);
            coupling.slave_surfaces = compiled_region(_data->surface_sets, coupling.slave_surfaces);
            _data->couplings.push_back(std::move(coupling));
        }

        for (auto rigid_body : source->rigid_bodies) {
            rigid_body.ref_node   = compiled_node(rigid_body.ref_node);
            rigid_body.ref_region = compiled_region(_data->node_sets, rigid_body.ref_region);
            rigid_body.elements   = compiled_region(_data->elem_sets, rigid_body.elements);
            rigid_body.tie_nodes  = compiled_region(_data->node_sets, rigid_body.tie_nodes);
            rigid_body.pin_nodes  = compiled_region(_data->node_sets, rigid_body.pin_nodes);
            _data->rigid_bodies.push_back(std::move(rigid_body));
        }

        // Geometry adjustment of copied ties is performed by the reader once the
        // initial configuration is final.
        for (auto tie : source->ties) {
            tie.master_surfaces = compiled_region(_data->surface_sets, tie.master_surfaces);
            tie.master_lines    = compiled_region(_data->line_sets   , tie.master_lines);
            tie.slave_nodes     = compiled_region(_data->node_sets   , tie.slave_nodes);
            tie.slave_surfaces  = compiled_region(_data->surface_sets, tie.slave_surfaces);
            _data->ties.push_back(std::move(tie));
        }

        for (auto contact : source->contacts) {
            contact.master_surfaces = compiled_region(_data->surface_sets, contact.master_surfaces);
            contact.slave_surfaces  = compiled_region(_data->surface_sets, contact.slave_surfaces);
            _data->contacts.push_back(std::move(contact));
        }

        for (auto rbm : source->rbms) {
            rbm.element_region = compiled_region(_data->elem_sets, rbm.element_region);
            _data->rbms.push_back(std::move(rbm));
        }

        // Equation DOFs keep their global directions, only the nodes are remapped
        for (auto equation : source->equations) {
            for (auto& entry : equation.entries) {
                entry.node_id = compiled_node(entry.node_id);
            }
            _data->equations.push_back(std::move(equation));
        }
    }

    // Section assignment and element-point enumeration operate exclusively on
    // the now dense assembly. Once these are complete, semantic topology is
    // frozen permanently and all later operations use global identifiers.
    assign_sections();
    _data->initialize_element_enumeration();
    _data->compiled = true;
}

/**
 * Resolves a compact node reference into dense assembly node space.
 *
 * The reference may contain a local identifier for the implicit identity
 * Instance or an explicitly qualified `INSTANCE.ID` value.
 *
 * @param reference Local `ID` or qualified `INSTANCE.ID` node reference.
 * @return Dense assembly node identifier.
 */
Index Model::compiled_node_id(const std::string& reference) const {
    const auto [instance, id] = split_compiled_reference(reference, "node");
    return compiled_node_id(id, instance);
}

/**
 * Resolves a compact element reference into dense assembly element space.
 *
 * The reference may contain a local identifier for the implicit identity
 * Instance or an explicitly qualified `INSTANCE.ID` value.
 *
 * @param reference Local `ID` or qualified `INSTANCE.ID` element reference.
 * @return Dense assembly element identifier.
 */
Index Model::compiled_element_id(const std::string& reference) const {
    const auto [instance, id] = split_compiled_reference(reference, "element");
    return compiled_element_id(id, instance);
}

/**
 * Resolves a compact surface reference into dense assembly surface space.
 *
 * The reference may contain a local identifier for the implicit identity
 * Instance or an explicitly qualified `INSTANCE.ID` value.
 *
 * @param reference Local `ID` or qualified `INSTANCE.ID` surface reference.
 * @return Dense assembly surface identifier.
 */
Index Model::compiled_surface_id(const std::string& reference) const {
    const auto [instance, id] = split_compiled_reference(reference, "surface");
    return compiled_surface_id(id, instance);
}

/**
 * Resolves a scalar line reference in an Instance namespace.
 *
 * @param reference Local `ID` or qualified `INSTANCE.ID` line reference.
 * @return Dense compiled line identifier.
 */
Index Model::compiled_line_id(const std::string& reference) const {
    const auto [instance, id] = split_compiled_reference(reference, "line");
    return compiled_line_id(id, instance);
}

/**
 * Resolves a node reference on the current side of the compile boundary.
 *
 * Part-scoped commands run before compilation and keep part-local node ids, so
 * the reference must be an unqualified integer there. Assembly-scoped commands
 * run afterwards and resolve `ID` or `INSTANCE.ID` through the Instance maps.
 *
 * @param reference Node reference as given in the input.
 * @return Part-local id before compile(), dense compiled id afterwards.
 */
Index Model::resolve_node_id(const std::string& reference) const {
    logging::error(_data != nullptr,
        "Model: model data is not initialized");

    if (_data->compiled) {
        return compiled_node_id(reference);
    }

    logging::error(reference.find('.') == std::string::npos,
        "Model: node reference '", reference, "' cannot be instance-qualified inside a part");
    return split_compiled_reference(reference, "node").second;
}

/**
 * Resolves one part-local node identifier through a compiled Instance map.
 *
 * @param id Node identifier in the referenced Part namespace.
 * @param instance Instance name, or the implicit identity instance by default.
 * @return Dense assembly node identifier.
 */
Index Model::compiled_node_id(Index id, const std::string& instance) const {
    // Validate the compile state, instance namespace and local node reference
    logging::error(_data != nullptr && _data->compiled,
        "Model: node ids cannot be resolved before compile()");

    const auto compiled_instance = _data->instances.get(instance);
    logging::error(compiled_instance != nullptr,
        "Model: compiled instance ", instance, " does not exist");

    const auto it = compiled_instance->node_ids.find(id);
    logging::error(it != compiled_instance->node_ids.end(),
        "Model: node ", id, " does not exist in instance ", instance);
    return it->second;
}

/**
 * Resolves one part-local element identifier through a compiled Instance map.
 *
 * @param id Element identifier in the referenced Part namespace.
 * @param instance Instance name, or the implicit identity instance by default.
 * @return Dense assembly element identifier.
 */
Index Model::compiled_element_id(Index id, const std::string& instance) const {
    // Validate the compile state, instance namespace and local element reference
    logging::error(_data != nullptr && _data->compiled,
        "Model: element ids cannot be resolved before compile()");

    const auto compiled_instance = _data->instances.get(instance);
    logging::error(compiled_instance != nullptr,
        "Model: compiled instance ", instance, " does not exist");

    const auto it = compiled_instance->element_ids.find(id);
    logging::error(it != compiled_instance->element_ids.end(),
        "Model: element ", id, " does not exist in instance ", instance);
    return it->second;
}

/**
 * Resolves one part-local surface identifier through a compiled Instance map.
 *
 * @param id Surface identifier in the referenced Part namespace.
 * @param instance Instance name, or the implicit identity instance by default.
 * @return Dense assembly surface identifier.
 */
Index Model::compiled_surface_id(Index id, const std::string& instance) const {
    // Validate the compile state, instance namespace and local surface reference
    logging::error(_data != nullptr && _data->compiled,
        "Model: surface ids cannot be resolved before compile()");

    const auto compiled_instance = _data->instances.get(instance);
    logging::error(compiled_instance != nullptr,
        "Model: compiled instance ", instance, " does not exist");

    const auto it = compiled_instance->surface_ids.find(id);
    logging::error(it != compiled_instance->surface_ids.end(),
        "Model: surface ", id, " does not exist in instance ", instance);
    return it->second;
}

/**
 * Resolves one part-local line identifier through a compiled Instance map.
 *
 * @param id Line identifier in the referenced Part namespace.
 * @param instance Instance name, or the implicit identity instance by default.
 * @return Dense assembly line identifier.
 */
Index Model::compiled_line_id(Index id, const std::string& instance) const {
    // Validate the compile state, instance namespace and local line reference
    logging::error(_data != nullptr && _data->compiled,
        "Model: line ids cannot be resolved before compile()");

    const auto compiled_instance = _data->instances.get(instance);
    logging::error(compiled_instance != nullptr,
        "Model: compiled instance ", instance, " does not exist");

    const auto it = compiled_instance->line_ids.find(id);
    logging::error(it != compiled_instance->line_ids.end(),
        "Model: line ", id, " does not exist in instance ", instance);
    return it->second;
}

/**
 * @brief Resolves a compiled node set or scalar node reference into a node region.
 *
 * Named node sets are returned directly. Scalar references may be given as an
 * unqualified node ID or as `INSTANCE.ID` and are mapped to the corresponding
 * compiled assembly identifier before a private single-node region is created.
 *
 * @param reference Named node set or scalar node reference.
 * @return Compiled node region represented by the reference.
 */
NodeRegion::Ptr Model::resolve_node_region(const std::string& reference) const {
    // Reuse an existing compiled node set whenever the reference names one.
    if (_data->node_sets.has(reference)) {
        return _data->node_sets.get(reference);
    }

    // Scalar node references are mapped into compiled assembly space and
    // represented by a private single-node region.
    auto region = std::make_shared<NodeRegion>("INTERNAL");
    region->add(compiled_node_id(reference));

    return region;
}

/**
 * @brief Resolves a compiled element set or scalar element reference into an element region.
 *
 * Named element sets are returned directly. Scalar references may be given as
 * an unqualified element ID or as `INSTANCE.ID` and are mapped to the corresponding
 * compiled assembly identifier before a private single-element region is created.
 *
 * @param reference Named element set or scalar element reference.
 * @return Compiled element region represented by the reference.
 */
ElementRegion::Ptr Model::resolve_element_region(const std::string& reference) const {
    // Reuse an existing compiled element set whenever the reference names one.
    if (_data->elem_sets.has(reference)) {
        return _data->elem_sets.get(reference);
    }

    // Scalar element references are mapped into compiled assembly space and
    // represented by a private single-element region.
    auto region = std::make_shared<ElementRegion>("INTERNAL");
    region->add(compiled_element_id(reference));

    return region;
}

/**
 * @brief Resolves a compiled surface set or scalar surface reference into a region.
 *
 * Named surface sets are returned directly. Scalar references may be given as
 * an unqualified surface ID or as INSTANCE.ID and are mapped to the compiled
 * assembly identifier before a private single-surface region is created.
 * Reference validation and instance mapping remain in compiled_surface_id().
 *
 * @param reference Named surface set or scalar surface reference.
 * @return Compiled surface region represented by the reference.
 */
SurfaceRegion::Ptr Model::resolve_surface_region(const std::string& reference) const {
    // Preserve shared named regions instead of copying their surface IDs.
    if (_data->surface_sets.has(reference)) {
        return _data->surface_sets.get(reference);
    }

    // Map scalar references into assembly space before creating a private region.
    auto region = std::make_shared<SurfaceRegion>("INTERNAL");
    region->add(compiled_surface_id(reference));
    return region;
}

} // namespace fem::model
