/**
 * @file model_compile.cpp
 * @brief Compiles part/instance topology into dense solver storage.
 *
 * `Model::compile()` is intentionally narrow: it freezes only the semantic
 * part/instance topology that must be flattened into the assembly. Nodes,
 * elements, surfaces and lines are copied from every instance, transformed and
 * rewired to dense global identifiers. Instance-local sets and section regions
 * are materialized in the same global identifier space.
 *
 * Definition objects such as materials, profiles and coordinate systems are not
 * part of that one-way topology transition. They may be registered before or
 * after compilation. Sections may likewise be added after compilation when they
 * already reference a compiled element set; pre-compile sections remain part
 * definitions and are copied once per instance during this routine.
 *
 * Section orientations require the same rigid placement as their owning part.
 * A compiled section therefore receives an independent transformed coordinate
 * system: orientation axes rotate with the instance and spatial origins also
 * receive the instance translation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "model.h"

#include "../section/section_beam.h"
#include "../section/section_shell_abd.h"
#include "../section/section_shell_integrated.h"
#include "../section/section_solid.h"
#include "../section/section_truss.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace fem::model {

Model::Model()
    : _data(std::make_shared<ModelData>()) {
    // Unqualified deck entities live in the default part. The corresponding
    // default instance is its identity embedding into the final assembly.
    const auto default_part = _data->parts.activate(DEFAULT_PART_NAME);
    _data->instances.activate(DEFAULT_INSTANCE_NAME, default_part);
}

void Model::add_part(const std::string& name) {
    // Parts participate directly in flattening and may therefore only be
    // created before the one allowed compile pass.
    logging::error(_data != nullptr && !_data->compiled,
        "Model: parts cannot be added after compile()");
    logging::error(!name.empty(),
        "Model: part name cannot be empty");
    logging::error(!_data->parts.has(name),
        "Model: part ", name, " is already defined");

    _data->parts.activate(name);
}

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

    _data->instances.add(std::make_shared<Instance>(
        name,
        _data->parts.get(part_name),
        translation,
        rotation
    ));
}

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

    _data->coordinate_systems.add(std::move(coordinate_system));
}

void Model::add_material(material::Material::Ptr material) {
    // Materials do not affect part/instance identifier flattening. They may be
    // introduced whenever a later section or analysis object needs them.
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(material != nullptr,
        "Model: cannot add a null material");
    logging::error(!_data->materials.has(material->name),
        "Model: material ", material->name, " is already defined");

    _data->materials.add(std::move(material));
}

void Model::add_profile(Profile::Ptr profile) {
    // Profiles behave like materials: they are shared definitions and are not
    // part of the topology that compile() freezes.
    logging::error(_data != nullptr,
        "Model: model data is not initialized");
    logging::error(profile != nullptr,
        "Model: cannot add a null profile");
    logging::error(!_data->profiles.has(profile->name),
        "Model: profile ", profile->name, " is already defined");

    _data->profiles.add(std::move(profile));
}

void Model::add_section(Section::Ptr section) {
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
        logging::error(_data->elem_sets.has(section->region_->name)
                    && _data->elem_sets.get(section->region_->name) == section->region_,
            "Model: post-compile section region ", section->region_->name,
            " is not a compiled element set");

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
    logging::error(active->elem_sets.has(section->region_->name)
                && active->elem_sets.get(section->region_->name) == section->region_,
        "Model: section region ", section->region_->name,
        " does not belong to active part ", active->name);

    active->sections.push_back(std::move(section));
}

void Model::compile() {
    // Compiling twice would duplicate/renumber the assembly and invalidate every
    // global reference created after the first pass. Recompilation is therefore
    // explicitly unsupported instead of attempting cache invalidation.
    logging::error(_data != nullptr && !_data->compiled,
        "Model: compile() may only be called once");

    std::vector<std::string> instance_names;
    instance_names.reserve(_data->instances._data.size());

    Index total_nodes    = 0;
    Index total_elements = 0;
    Index total_surfaces = 0;
    Index total_lines    = 0;

    // Only proper rigid rotations are accepted. Reflections would change element
    // orientation and surface handedness and are not valid instance placements.
    const Precision rotation_tolerance =
        Precision(100) * std::sqrt(std::numeric_limits<Precision>::epsilon());

    for (const auto& [name, instance] : _data->instances) {
        logging::error(instance != nullptr && instance->part != nullptr,
            "Model: instance ", name, " does not reference a part");
        logging::error(instance->translation.allFinite(),
            "Model: instance ", name, " translation contains invalid values");
        logging::error(instance->rotation.allFinite(),
            "Model: instance ", name, " rotation contains invalid values");
        logging::error((instance->rotation.transpose() * instance->rotation - Mat3::Identity()).norm()
                    <= rotation_tolerance,
            "Model: instance ", name, " rotation is not orthonormal");
        logging::error(std::abs(instance->rotation.determinant() - Precision(1)) <= rotation_tolerance,
            "Model: instance ", name, " rotation must have determinant +1");

        // The local-to-global maps are output of compile(). Clearing them here
        // makes the single pass deterministic even for a freshly reused object,
        // while the compiled guard still forbids an actual second invocation.
        instance->node_ids.clear();
        instance->element_ids.clear();
        instance->surface_ids.clear();
        instance->line_ids.clear();

        total_nodes    += static_cast<Index>(instance->part->nodes.size());
        total_elements += static_cast<Index>(instance->part->elements.size());
        total_surfaces += static_cast<Index>(instance->part->surfaces.size());
        total_lines    += static_cast<Index>(instance->part->lines.size());
        instance_names.push_back(name);
    }
    std::sort(instance_names.begin(), instance_names.end());

    // Dense solver arrays are allocated exactly once from the semantic assembly
    // size. Sparse/local identifiers never escape into these arrays.
    _data->elements.assign(static_cast<std::size_t>(total_elements), nullptr);
    _data->surfaces.assign(static_cast<std::size_t>(total_surfaces), nullptr);
    _data->lines.assign(static_cast<std::size_t>(total_lines), nullptr);
    _data->sections.clear();

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
    _data->positions = std::make_shared<Field>(
        "POSITION", FieldDomain::NODE, total_nodes, 6);
    _data->positions_reference = std::make_shared<Field>(
        "POSITION_REFERENCE", FieldDomain::NODE, total_nodes, 6);
    _data->positions->set_zero();
    _data->positions_reference->set_zero();
    _data->fields.emplace(_data->positions->name, _data->positions);
    _data->fields.emplace(_data->positions_reference->name, _data->positions_reference);

    ID next_node    = 0;
    ID next_element = 0;
    ID next_surface = 0;
    ID next_line    = 0;

    for (const std::string& instance_name : instance_names) {
        const auto instance   = _data->instances.get(instance_name);
        const auto source     = instance->part;
        const bool is_default = instance_name == DEFAULT_INSTANCE_NAME;

        // The default instance owns the unqualified assembly namespace. Explicit
        // instances prefix inherited part sets with "INSTANCE.".
        auto qualified_name = [&](const std::string& local_name) {
            return is_default ? local_name : instance_name + "." + local_name;
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

        // Sort semantic identifiers before assigning dense ids so compiled output
        // is deterministic despite unordered_map storage in Part.
        std::vector<ID> node_ids;
        node_ids.reserve(source->nodes.size());
        for (const auto& [id, position] : source->nodes) {
            (void) position;
            node_ids.push_back(id);
        }
        std::sort(node_ids.begin(), node_ids.end());

        for (const ID local_id : node_ids) {
            const ID global_id = next_node++;
            const Vec3 position =
                instance->rotation * source->nodes.at(local_id) + instance->translation;

            node_map.emplace(local_id, global_id);
            for (Dim d = 0; d < 3; ++d) {
                (*_data->positions)(static_cast<Index>(global_id), d)           = position(d);
                (*_data->positions_reference)(static_cast<Index>(global_id), d) = position(d);
            }
            _data->node_sets.all()->add(global_id);
        }

        std::vector<ID> element_ids;
        element_ids.reserve(source->elements.size());
        for (const auto& [id, element] : source->elements) {
            (void) element;
            element_ids.push_back(id);
        }
        std::sort(element_ids.begin(), element_ids.end());

        for (const ID local_id : element_ids) {
            const auto source_element = source->elements.at(local_id);
            const ID global_id        = next_element++;
            auto element              = source_element->copy();

            // copy() preserves the concrete element type and topology-specific
            // configuration. Assembly-owned ids, offsets, section assignment and
            // ModelData binding are always rebuilt explicitly here.
            element->elem_id           = global_id;
            element->elem_nodal_offset = 0;
            element->elem_ip_offset    = 0;
            element->elem_mp_offset    = 0;
            element->_section          = nullptr;
            element->_model_data       = _data.get();

            for (Index local = 0; local < static_cast<Index>(element->n_nodes()); ++local) {
                const ID local_node = source_element->nodes()[local];
                const auto node_it  = node_map.find(local_node);
                logging::error(node_it != node_map.end(),
                    "Model: element ", local_id, " in part ", source->name,
                    " references undefined node ", local_node);
                element->nodes()[local] = node_it->second;
            }

            _data->elements[static_cast<std::size_t>(global_id)] = std::move(element);
            _data->elem_sets.all()->add(global_id);
            element_map.emplace(local_id, global_id);
        }

        // Materialize all inherited part-local node and element sets. The
        // default *ALL sets are already global parents and must not be replaced by
        // a second same-named default-part collection.
        for (const auto& [local_name, local_set] : source->node_sets) {
            if (is_default && local_name == SET_NODE_ALL) continue;

            auto target = _data->node_sets.activate(qualified_name(local_name));
            for (const ID local_id : *local_set) {
                const auto it = node_map.find(local_id);
                logging::error(it != node_map.end(),
                    "Model: node set ", local_name, " in part ", source->name,
                    " references undefined node ", local_id);
                target->add(it->second);
            }
        }

        for (const auto& [local_name, local_set] : source->elem_sets) {
            if (is_default && local_name == SET_ELEM_ALL) continue;

            auto target = _data->elem_sets.activate(qualified_name(local_name));
            for (const ID local_id : *local_set) {
                const auto it = element_map.find(local_id);
                logging::error(it != element_map.end(),
                    "Model: element set ", local_name, " in part ", source->name,
                    " references undefined element ", local_id);
                target->add(it->second);
            }
        }

        std::vector<ID> surface_ids;
        surface_ids.reserve(source->surfaces.size());
        for (const auto& [id, surface] : source->surfaces) {
            (void) surface;
            surface_ids.push_back(id);
        }
        std::sort(surface_ids.begin(), surface_ids.end());

        for (const ID local_id : surface_ids) {
            const auto source_surface = source->surfaces.at(local_id);
            const ID global_id        = next_surface++;
            auto surface              = source_surface->copy();

            for (Index local = 0; local < surface->n_nodes; ++local) {
                const ID local_node = source_surface->nodes()[local];
                const auto node_it  = node_map.find(local_node);
                logging::error(node_it != node_map.end(),
                    "Model: surface ", local_id, " in part ", source->name,
                    " references undefined node ", local_node);
                surface->nodes()[local] = node_it->second;
            }

            _data->surfaces[static_cast<std::size_t>(global_id)] = std::move(surface);
            _data->surface_sets.all()->add(global_id);
            surface_map.emplace(local_id, global_id);
        }

        for (const auto& [local_name, local_set] : source->surface_sets) {
            if (is_default && local_name == SET_SURF_ALL) continue;

            auto target = _data->surface_sets.activate(qualified_name(local_name));
            for (const ID local_id : *local_set) {
                const auto it = surface_map.find(local_id);
                logging::error(it != surface_map.end(),
                    "Model: surface set ", local_name, " in part ", source->name,
                    " references undefined surface ", local_id);
                target->add(it->second);
            }
        }

        std::vector<ID> line_ids;
        line_ids.reserve(source->lines.size());
        for (const auto& [id, line] : source->lines) {
            (void) line;
            line_ids.push_back(id);
        }
        std::sort(line_ids.begin(), line_ids.end());

        for (const ID local_id : line_ids) {
            const auto source_line = source->lines.at(local_id);
            const ID global_id     = next_line++;
            auto line              = source_line->copy();

            for (Index local = 0; local < line->n_nodes; ++local) {
                const ID local_node = source_line->nodes()[local];
                const auto node_it  = node_map.find(local_node);
                logging::error(node_it != node_map.end(),
                    "Model: line ", local_id, " in part ", source->name,
                    " references undefined node ", local_node);
                line->nodes()[local] = node_it->second;
            }

            _data->lines[static_cast<std::size_t>(global_id)] = std::move(line);
            _data->line_sets.all()->add(global_id);
            line_map.emplace(local_id, global_id);
        }

        for (const auto& [local_name, local_set] : source->line_sets) {
            if (is_default && local_name == SET_LINE_ALL) continue;

            auto target = _data->line_sets.activate(qualified_name(local_name));
            for (const ID local_id : *local_set) {
                const auto it = line_map.find(local_id);
                logging::error(it != line_map.end(),
                    "Model: line set ", local_name, " in part ", source->name,
                    " references undefined line ", local_id);
                target->add(it->second);
            }
        }

        // Sections are definitions on the source part but assignments in the
        // compiled assembly. Their regions and spatial orientations therefore
        // need instance-specific copies.
        for (const Section::Ptr& source_section : source->sections) {
            auto region = std::make_shared<ElementRegion>(
                qualified_name(source_section->region_->name));

            for (const ID local_id : *source_section->region_) {
                const auto it = element_map.find(local_id);
                logging::error(it != element_map.end(),
                    "Model: section region ", source_section->region_->name,
                    " in part ", source->name,
                    " references undefined element ", local_id);
                region->add(it->second);
            }

            Section::Ptr section = nullptr;
            if (auto* solid = source_section->as<SolidSection>()) {
                auto copy = std::make_shared<SolidSection>();
                copy->material_    = solid->material_;
                copy->region_      = region;
                copy->orientation_ = transformed_orientation(solid->orientation_);
                section = std::move(copy);
            } else if (auto* beam = source_section->as<BeamSection>()) {
                auto copy = std::make_shared<BeamSection>();
                copy->material_  = beam->material_;
                copy->region_    = region;
                copy->direction_ = instance->rotation * beam->direction_;
                copy->profile_   = beam->profile_;
                section = std::move(copy);
            } else if (auto* truss = source_section->as<TrussSection>()) {
                section = std::make_shared<TrussSection>(
                    truss->material_, region, truss->area_);
            } else if (auto* shell = source_section->as<IntegratedShellSection>()) {
                section = std::make_shared<IntegratedShellSection>(
                    shell->material_,
                    region,
                    shell->thickness_,
                    transformed_orientation(shell->orientation_),
                    shell->csys_axis_
                );
            } else if (auto* shell = source_section->as<ABDShellSection>()) {
                section = std::make_shared<ABDShellSection>(
                    shell->material_,
                    region,
                    shell->thickness_,
                    shell->abd_,
                    shell->shear_,
                    transformed_orientation(shell->orientation_),
                    shell->csys_axis_
                );
            }

            logging::error(section != nullptr,
                "Model: unsupported section type in part ", source->name);
            _data->sections.push_back(std::move(section));
        }
    }

    // Section assignment and element-point enumeration operate exclusively on
    // the now dense assembly. Once these are complete, semantic topology is
    // frozen permanently and all later operations use global identifiers.
    assign_sections();
    _data->initialize_element_enumeration();
    _data->compiled = true;
}

ID Model::compiled_node_id(ID id, const std::string& instance) const {
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

ID Model::compiled_element_id(ID id, const std::string& instance) const {
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

ID Model::compiled_surface_id(ID id, const std::string& instance) const {
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

ID Model::compiled_line_id(ID id, const std::string& instance) const {
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

} // namespace fem::model
