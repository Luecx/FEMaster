#include "model.h"
#include "element/element_structural.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace fem::model {

void Model::add_pretension_interface_section(const std::string& name,
                                             const std::string& surface_set_a,
                                             const std::string& surface_set_b) {
    logging::error(!name.empty(), "PRETENSION SECTION: name must not be empty");
    logging::error(surface_set_a != surface_set_b,
        "PRETENSION SECTION: two separate surface sets are required");
    logging::error(_data->surface_sets.has(surface_set_a) && _data->surface_sets.has(surface_set_b),
        "PRETENSION SECTION: both surface sets must exist");
    logging::error(std::none_of(_data->pretension_sections.begin(), _data->pretension_sections.end(),
        [&](const auto& section) { return section->name == name; }),
        "PRETENSION SECTION: duplicate section name ", name);
    auto section = std::make_shared<pretension::PretensionSection>();
    section->name = name;
    section->interface_surface_set_a = surface_set_a;
    section->interface_surface_set_b = surface_set_b;
    _data->pretension_sections.push_back(std::move(section));
}

namespace {
std::set<ID> surface_nodes(ModelData& data, const std::string& name) {
    std::set<ID> nodes;
    const auto region = data.surface_sets.get(name);
    logging::error(region && region->size() != 0, "PRETENSION SECTION: empty surface set ", name);
    for (const ID id : *region) {
        logging::error(id >= 0 && id < static_cast<ID>(data.surfaces.size()) && data.surfaces[id],
            "PRETENSION SECTION: invalid surface in ", name);
        logging::error(id < static_cast<ID>(data.surface_element_ids.size()) &&
            data.surface_element_ids[id] >= 0 &&
            data.surface_element_ids[id] < static_cast<ID>(data.elements.size()),
            "PRETENSION SECTION: surface has no owning element");
        const auto* element = data.elements[data.surface_element_ids[id]]->as<StructuralElement>();
        logging::error(element && element->is_solid(),
            "PRETENSION SECTION: only solid element faces are supported");
        const auto& surface = data.surfaces[id];
        for (Index i = 0; i < surface->n_nodes; ++i) nodes.insert(surface->nodes()[i]);
    }
    return nodes;
}
}

void Model::prepare_pretension_sections() {
    const auto& reference = *_data->positions_reference;
    for (auto& section : _data->pretension_sections) {
        const auto a = surface_nodes(*_data, section->interface_surface_set_a);
        const auto b = surface_nodes(*_data, section->interface_surface_set_b);
        logging::error(!a.empty() && a.size() == b.size(),
            "PRETENSION SECTION '", section->name, "': nonconforming surface node counts");
        for (ID node : a) logging::error(!b.count(node),
            "PRETENSION SECTION '", section->name, "': surfaces share topology; split the mesh in the preprocessor");
        for (const auto& previous : _data->pretension_sections) {
            if (previous == section || !previous->prepared) continue;
            for (const auto& pair : previous->interface_pairs)
                logging::error(!a.count(pair.side_a) && !a.count(pair.side_b) &&
                    !b.count(pair.side_a) && !b.count(pair.side_b),
                    "PRETENSION SECTION: overlapping sections are unsupported");
        }

        Vec3 normal = Vec3::Zero();
        for (ID id : *_data->surface_sets.get(section->interface_surface_set_a)) {
            const auto& surface = _data->surfaces[id];
            normal += surface->area(reference) * surface->normal(reference, Vec2::Zero());
        }
        logging::error(normal.allFinite() && normal.norm() > Precision(0),
            "PRETENSION SECTION: cannot determine cut normal");
        section->axis_direction = normal.normalized();
        for (ID id : *_data->surface_sets.get(section->interface_surface_set_a))
            logging::error(section->axis_direction.dot(
                _data->surfaces[id]->normal(reference, Vec2::Zero())) > Precision(0.999999),
                "PRETENSION SECTION: inconsistent cut normals");

        Vec3 center = Vec3::Zero();
        for (ID node : a) center += reference.row_vec3(node);
        center /= static_cast<Precision>(a.size());
        Precision scale = Precision(1);
        for (ID node : a) scale = std::max(scale, (reference.row_vec3(node) - center).norm());
        const Precision tolerance = Precision(1e-8) * scale;
        for (ID node : a)
            logging::error(std::abs(section->axis_direction.dot(reference.row_vec3(node) - center)) <= tolerance,
                "PRETENSION SECTION: cut must be planar");

        // ponytail: quadratic matching; use a spatial index if large cuts make preparation slow.
        std::set<ID> used_b;
        std::vector<pretension::InterfacePair> pairs;
        for (ID node_a : a) {
            ID match = -1;
            for (ID node_b : b) {
                if ((reference.row_vec3(node_a) - reference.row_vec3(node_b)).norm() > tolerance) continue;
                logging::error(match < 0,
                    "PRETENSION SECTION '", section->name, "': ambiguous match for A node ", node_a);
                match = node_b;
            }
            logging::error(match >= 0,
                "PRETENSION SECTION '", section->name, "': missing match for A node ", node_a);
            logging::error(used_b.insert(match).second,
                "PRETENSION SECTION '", section->name, "': B node has multiple matches ", match);
            pairs.push_back({node_a, match});
        }
        logging::error(used_b.size() == b.size(), "PRETENSION SECTION: unmatched B nodes");

        // A node bijection alone does not establish a conforming face mesh.
        std::map<ID, ID> matched_nodes;
        for (const auto& pair : pairs) matched_nodes.emplace(pair.side_a, pair.side_b);
        const auto face_nodes = [&](ID id, bool map_to_b) {
            std::vector<ID> nodes;
            const auto& face = _data->surfaces[id];
            for (Index i = 0; i < face->n_nodes; ++i) {
                const ID node = face->nodes()[i];
                nodes.push_back(map_to_b ? matched_nodes.at(node) : node);
            }
            std::sort(nodes.begin(), nodes.end());
            return nodes;
        };
        std::multiset<std::vector<ID>> faces_b;
        for (ID id : *_data->surface_sets.get(section->interface_surface_set_b))
            faces_b.insert(face_nodes(id, false));
        for (ID id : *_data->surface_sets.get(section->interface_surface_set_a)) {
            const auto match = faces_b.find(face_nodes(id, true));
            logging::error(match != faces_b.end(),
                "PRETENSION SECTION '", section->name, "': nonconforming face topology");
            faces_b.erase(match);
        }
        logging::error(faces_b.empty(),
            "PRETENSION SECTION: unmatched faces on SURFACE_B");
        section->interface_pairs = std::move(pairs);
        section->prepared = true;
    }
}

} // namespace fem::model
