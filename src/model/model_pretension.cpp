#include "model.h"
#include "element/element_structural.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>

namespace fem::model {

void Model::add_pretension_interface_section(const std::string& name,
                                             const std::string& surface_set_a,
                                             const std::string& surface_set_b) {
    logging::error(!name.empty(), "PRETENSION SECTION: name must not be empty");
    logging::error(_data->surface_sets.has(surface_set_a),
        "PRETENSION SECTION: surface set ", surface_set_a, " does not exist");
    logging::error(_data->surface_sets.has(surface_set_b),
        "PRETENSION SECTION: surface set ", surface_set_b, " does not exist");
    const bool duplicate = std::any_of(
        _data->pretension_sections.begin(), _data->pretension_sections.end(),
        [&](const auto& section) { return section && section->name == name; });
    logging::error(!duplicate,
        "PRETENSION SECTION: duplicate section name ", name);

    auto section = std::make_shared<pretension::PretensionSection>();
    section->name = name;
    section->interface_surface_set_a = surface_set_a;
    section->interface_surface_set_b = surface_set_b;
    _data->pretension_sections.push_back(std::move(section));
}

namespace {

std::set<ID> surface_nodes(ModelData& data, const std::string& set_name) {
    std::set<ID> nodes;
    for (const ID surface_id : *data.surface_sets.get(set_name)) {
        logging::error(surface_id >= 0 && surface_id < static_cast<ID>(data.surfaces.size()),
            "PRETENSION SECTION: invalid surface in ", set_name);
        const auto& surface = data.surfaces[static_cast<std::size_t>(surface_id)];
        logging::error(surface != nullptr, "PRETENSION SECTION: null surface in ", set_name);
        for (Index i = 0; i < surface->n_nodes; ++i) {
            nodes.insert(surface->nodes()[i]);
        }
    }
    return nodes;
}

std::set<ID> surface_owners(ModelData& data, const std::string& set_name) {
    std::set<ID> owners;
    for (const ID surface_id : *data.surface_sets.get(set_name)) {
        if (surface_id < 0 ||
            surface_id >= static_cast<ID>(data.surface_element_ids.size())) continue;
        const ID owner = data.surface_element_ids[static_cast<std::size_t>(surface_id)];
        if (owner >= 0) owners.insert(owner);
    }
    return owners;
}

// Geometry is checked on the original mesh, before any cut nodes are appended.
void validate_interface_geometry(ModelData& data,
                                 pretension::PretensionSection& section,
                                 const std::set<ID>& nodes_a,
                                 const Vec3& center_a, const Vec3& center_b,
                                 Precision tolerance) {
    Vec3 axis = center_b - center_a;
    if (axis.norm() <= tolerance) {
        SurfacePtr reference;
        for (const ID surface_id : *data.surface_sets.get(
                 section.interface_surface_set_a)) {
            if (surface_id >= 0 &&
                surface_id < static_cast<ID>(data.surfaces.size())) {
                reference = data.surfaces[static_cast<std::size_t>(surface_id)];
                if (reference) break;
            }
        }
        logging::error(reference != nullptr,
            "PRETENSION SECTION '", section.name,
            "': SURFACE_A contains no valid surface");
        axis = reference->normal(*data.positions, Vec2::Zero());
    }
    logging::error(axis.allFinite() && axis.norm() > Precision(0),
        "PRETENSION SECTION '", section.name,
        "': cannot determine interface axis");
    section.axis_direction = axis.normalized();
    for (const ID node : nodes_a) {
        logging::error(std::abs(section.axis_direction.dot(
            data.positions->row_vec3(static_cast<Index>(node)) - center_a)) <= tolerance,
            "PRETENSION SECTION '", section.name, "': interface must be planar");
    }
    if (section.nodal_control) {
        // Abaqus positive load acts on the underlying (rewired B) side along
        // the outward average surface normal. Native q uses the same B-A sign.
        Vec3 normal = Vec3::Zero();
        for (const ID id : *data.surface_sets.get(section.interface_surface_set_a)) {
            const auto& surface = data.surfaces[id];
            normal += surface->area(*data.positions) * surface->normal(*data.positions, Vec2::Zero());
        }
        logging::error(normal.allFinite() && normal.norm() > Precision(0),
            "PRE-TENSION SECTION: cannot determine average normal");
        normal.normalize();
        for (const ID id : *data.surface_sets.get(section.interface_surface_set_a))
            logging::error(normal.dot(data.surfaces[id]->normal(*data.positions, Vec2::Zero())) > Precision(0.999999),
                "PRE-TENSION SECTION: inconsistent surface orientation");
        section.axis_direction = section.normal_override.isZero()
            ? normal : section.normal_override.normalized();
    }

}

void prepare_surface_pair(ModelData& data,
                          pretension::PretensionSection& section,
                          const std::vector<std::vector<ID>>& node_elements) {
    for (const auto& set_name : {section.interface_surface_set_a, section.interface_surface_set_b}) {
        for (const ID surface : *data.surface_sets.get(set_name)) {
            logging::error(surface >= 0 && surface < static_cast<ID>(data.surface_element_ids.size()),
                "PRETENSION SECTION: surface has no owning solid element");
            const ID owner = data.surface_element_ids[surface];
            logging::error(owner >= 0 && owner < static_cast<ID>(data.elements.size()),
                "PRETENSION SECTION: invalid surface owner");
            const auto* structural = data.elements[owner]->as<StructuralElement>();
            logging::error(structural && structural->is_solid(),
                "PRETENSION SECTION: only solid element faces are supported");
        }
    }
    const auto nodes_a = surface_nodes(data, section.interface_surface_set_a);
    const auto nodes_b = surface_nodes(data, section.interface_surface_set_b);
    logging::error(!nodes_a.empty() && !nodes_b.empty(),
        "PRETENSION SECTION '", section.name,
        "': both interface surfaces must contain nodes");
    logging::error(nodes_a.size() == nodes_b.size(),
        "PRETENSION SECTION '", section.name,
        "': SURFACE_A and SURFACE_B have different node counts");

    Vec3 center_a = Vec3::Zero();
    Vec3 center_b = Vec3::Zero();
    Precision scale = 0;
    for (const ID node : nodes_a) {
        const Vec3 point = data.positions->row_vec3(static_cast<Index>(node));
        center_a += point;
    }
    for (const ID node : nodes_b) {
        const Vec3 point = data.positions->row_vec3(static_cast<Index>(node));
        center_b += point;
    }
    center_a /= static_cast<Precision>(nodes_a.size());
    center_b /= static_cast<Precision>(nodes_b.size());
    for (const ID node : nodes_a) {
        scale = std::max(scale,
            (data.positions->row_vec3(static_cast<Index>(node)) - center_a).norm());
    }
    const Precision tolerance = Precision(1e-8) * std::max(scale, Precision(1));
    for (const auto& existing : data.pretension_sections) {
        if (!existing || !existing->prepared) continue;
        for (const auto& pair : existing->interface_pairs) {
            logging::error(!nodes_a.count(pair.side_a) && !nodes_a.count(pair.side_b) &&
                          !nodes_b.count(pair.side_a) && !nodes_b.count(pair.side_b),
                "PRETENSION SECTION '", section.name, "': overlapping sections are unsupported");
        }
    }

    // Match coincident nodes; keep deterministic nearest-node/tie ordering.
    // ponytail: quadratic matching; add spatial search if larger cuts justify it.
    std::map<ID, ID> pairs;
    std::set<ID> unused_b = nodes_b;
    for (const ID node_a : nodes_a) {
        const Vec3 point_a = data.positions->row_vec3(static_cast<Index>(node_a));
        ID best = -1;
        Precision distance_best = std::numeric_limits<Precision>::max();
        for (const ID node_b : unused_b) {
            const Precision distance = (point_a - data.positions->row_vec3(
                static_cast<Index>(node_b))).norm();
            if (distance < distance_best) {
                distance_best = distance;
                best = node_b;
            }
        }
        logging::error(best >= 0 && distance_best <= tolerance,
            "PRETENSION SECTION '", section.name,
            "': interface nodes cannot be paired within tolerance ", tolerance);
        pairs.emplace(node_a, best);
        unused_b.erase(best);
    }
    logging::error(unused_b.empty(),
        "PRETENSION SECTION '", section.name,
        "': not all SURFACE_B nodes could be paired");

    validate_interface_geometry(data, section, nodes_a, center_a, center_b, tolerance);

    std::set<ID> shared_nodes;
    for (const auto& [node_a, node_b] : pairs) {
        if (node_a == node_b) shared_nodes.insert(node_a);
    }

    std::set<ID> side_b_elements;
    std::map<ID, ID> duplicated_b;
    if (!shared_nodes.empty()) {
        const auto seeds = surface_owners(data, section.interface_surface_set_b);
        logging::error(!seeds.empty(),
            "PRETENSION SECTION '", section.name,
            "': SURFACE_B has no owning solid elements");

        side_b_elements = seeds;
        std::deque<ID> queue(seeds.begin(), seeds.end());
        while (!queue.empty()) {
            const ID element_id = queue.front();
            queue.pop_front();
            const auto& element = data.elements[static_cast<std::size_t>(element_id)];
            for (const ID node : *element) {
                if (shared_nodes.count(node)) continue;
                for (const ID neighbor : node_elements[node]) {
                    if (side_b_elements.insert(neighbor).second) queue.push_back(neighbor);
                }
            }
        }
        for (const ID node : shared_nodes) {
            bool retained = false;
            for (const ID element : node_elements[node])
                retained |= !side_b_elements.count(element);
            logging::error(retained, "PRETENSION SECTION '", section.name,
                "': surface does not separate two element regions");
        }
        // Check every affected element before the first topology mutation.
        for (const ID node : shared_nodes)
            for (const ID id : node_elements[node]) {
                if (!side_b_elements.count(id)) continue;
                const auto* solid = data.elements[id]->as<StructuralElement>();
                logging::error(solid && solid->is_solid(),
                    "PRETENSION SECTION: splitting requires supported solid elements");
            }
    }

    // Allocate the cut copies and the native scalar control in one batch.
    std::vector<Vec3> new_positions;
    new_positions.reserve(shared_nodes.size() + (section.control_node < 0));
    for (const ID node : shared_nodes)
        new_positions.push_back(data.positions->row_vec3(node));
    if (section.control_node < 0) new_positions.push_back(Vec3::Zero());
    if (!new_positions.empty()) {
        const ID first = data.append_nodes(new_positions);
        ID next = first;
        for (const ID node : shared_nodes) duplicated_b.emplace(node, next++);
        if (section.control_node < 0) section.control_node = next;
    }

    // Rewire only incident elements, then let all owned boundary faces follow.
    if (!shared_nodes.empty()) {
        for (const auto& [old_node, new_node] : duplicated_b)
            for (const ID element_id : node_elements[old_node]) {
                if (!side_b_elements.count(element_id)) continue;
                logging::error(data.elements[element_id]->replace_node(old_node, new_node),
                    "PRETENSION SECTION: splitting requires supported solid elements");
            }
        // Every boundary on rewired elements must follow the new connectivity,
        // including surfaces later used by pressure, tie or contact definitions.
        for (ID surface_id = 0; surface_id < static_cast<ID>(data.surfaces.size()); ++surface_id) {
            if (!side_b_elements.count(data.surface_element_ids[surface_id])) continue;
            auto& surface = data.surfaces[static_cast<std::size_t>(surface_id)];
            if (!surface) continue;
            for (Index i = 0; i < surface->n_nodes; ++i) {
                const auto found = duplicated_b.find(surface->nodes()[i]);
                if (found != duplicated_b.end()) surface->nodes()[i] = found->second;
            }
        }
    }

    section.interface_pairs.clear();
    for (const auto& [node_a, original_b] : pairs) {
        const auto duplicated = duplicated_b.find(original_b);
        const ID node_b = duplicated == duplicated_b.end()
            ? original_b : duplicated->second;
        section.interface_pairs.push_back({node_a, node_b});
    }

    logging::info(true, "PRETENSIONSECTION '", section.name,
        "': paired ", section.interface_pairs.size(),
        " interface node(s), split ", shared_nodes.size(), " shared node(s)");
}

} // namespace

void Model::prepare_pretension_sections() {
    for (auto& section : _data->pretension_sections) {
        logging::error(section != nullptr,
            "PRETENSION SECTION: null section encountered");
        if (section->prepared) continue;
        // Rebuild per section: an earlier independent cut may have changed IDs.
        std::vector<std::vector<ID>> node_elements(_data->positions->rows);
        for (const auto& element : _data->elements) {
            if (!element) continue;
            for (const ID node : *element) node_elements[node].push_back(element->elem_id);
        }
        if (section->control_node >= 0) {
            logging::error(section->control_node < _data->positions->rows,
                "PRETENSION SECTION: invalid control node");
            for (const auto& other : _data->pretension_sections)
                logging::error(other == section || other->control_node != section->control_node,
                    "PRETENSION SECTION: control node used by multiple sections");
            logging::error(node_elements[section->control_node].empty(),
                "PRETENSION SECTION: control node must not be attached to an element");
        }
        prepare_surface_pair(*_data, *section, node_elements);
        section->prepared = true;
    }
}

void Model::set_pretension_load(const std::string& name,
                                pretension::Control control,
                                Precision value) {
    logging::error(std::isfinite(value), "PRETENSION: VALUE must be finite");
    for (auto& section : _data->pretension_sections) {
        if (!section || section->name != name) continue;
        logging::error(!section->nodal_control,
            "PRETENSION: native actions cannot be mixed with Abaqus control-node conditions");
        section->control = control;
        section->prescribed_value = value;
        section->state = pretension::State::Loading;
        section->has_solved_gap = false;
        return;
    }
    logging::error(false, "PRETENSION: section ", name, " does not exist");
}

void Model::lock_pretension_section(const std::string& name) {
    for (auto& section : _data->pretension_sections) {
        if (!section || section->name != name) continue;
        logging::error(!section->nodal_control,
            "PRETENSION: use BOUNDARY,FIXED for an Abaqus section");
        logging::error(section->state != pretension::State::Open,
            "PRETENSION LOCK: section ", name, " has not been loaded yet");
        if (section->state != pretension::State::Locked) {
            logging::error(section->has_solved_gap,
                "PRETENSION LOCK: section ", name, " has no solved gap to lock");
            section->locked_gap = section->last_solved_gap;
            section->control = pretension::Control::Displacement;
        }
        section->state = pretension::State::Locked;
        logging::info(true, "PRETENSION '", name,
            "': locked gap = ", section->locked_gap);
        return;
    }
    logging::error(false, "PRETENSION LOCK: section ", name, " does not exist");
}

} // namespace fem::model
