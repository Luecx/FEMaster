#include "element/element.h"
#include "element/element_structural.h"
#include "geometry/surface/surface.h"

inline void Model::set_node(ID id, Precision x, Precision y, Precision z) {
    logging::error(id < _data->max_nodes, "internal error; allocated less data than required. id=", id, " exceeds maximum limit");
    auto& node_coords = _data->node_data.get(NodeDataEntries::POSITION);

    logging::error(node_coords(id, 0) == 0 &&
                   node_coords(id, 1) == 0 &&
                   node_coords(id, 2) == 0, "node with id=", id, " seems to define non-zero values twice");

    node_coords(id, 0) = x;
    node_coords(id, 1) = y;
    node_coords(id, 2) = z;
    _data->node_sets.add(id);
}

template<typename T, typename... Args>
inline void Model::set_element(ID id, Args&&... args) {
    logging::error(id < _data->max_elems, "internal error; allocated less data than required. id=", id, " exceeds maximum limit");

    auto el = ElementPtr{new T{id, {args...}}};
    if (element_dims == 0) {
        element_dims = el->dimensions();
    }
    el->_model_data = _data;
    logging::error(element_dims == el->dimensions(), "different number of dimensions across elements of model; element with id=",
                   id, " has ", el->dimensions(), " while all other elements have ", element_dims, " dimensions");

    logging::error(_data->elements[id] == nullptr, "element with id=", id, " has already been defined");

    _data->elements[id] = std::move(el);
    _data->elem_sets.add(id);
}

template<typename T, typename... Args>
inline void Model::set_beam_element(ID id, ID orientation_node, Args&&... args) {
    logging::error(id < _data->max_elems, "internal error; allocated less data than required. id=", id, " exceeds maximum limit");

    auto el = ElementPtr{new T{id, {args...}, orientation_node}};
    if (element_dims == 0) {
        element_dims = el->dimensions();
    }
    el->_model_data = _data;
    logging::error(element_dims == el->dimensions(), "different number of dimensions across elements of model; element with id=",
                   id, " has ", el->dimensions(), " while all other elements have ", element_dims, " dimensions");

    logging::error(_data->elements[id] == nullptr, "element with id=", id, " has already been defined");

    _data->elements[id] = std::move(el);
    _data->elem_sets.add(id);
}

inline void Model::set_surface(ID id, ID element_id, ID surface_id) {
    logging::error(id < _data->max_surfaces, "internal error; allocated less data than required. id=", id, " exceeds maximum limit");

    auto& elptr = _data->elements[element_id];
    logging::error(elptr != nullptr, "element with id=", element_id, " has not been defined");
    auto surfptr = elptr->surface(surface_id);
    // If the element has no surface for the given side, consider 1D geometry case
    // (e.g., beams or trusses with 2-node connectivity) and create a line entry.
    if (surfptr == nullptr) {
        // Heuristic: elements with exactly 2 nodes represent 1D members here
        if (elptr->n_nodes() == 2) {
            // allow negative id = automatically assign id for lines too
            std::array<ID,2> ln{elptr->nodes()[0], elptr->nodes()[1]};

            // Deduplicate: reuse existing line id if same connectivity already present (either order)
            ID existing = -1;
            for (std::size_t i = 0; i < _data->lines.size(); ++i) {
                auto cur = _data->lines[i];
                if ((cur[0] == ln[0] && cur[1] == ln[1]) || (cur[0] == ln[1] && cur[1] == ln[0])) {
                    existing = static_cast<ID>(i);
                    break;
                }
            }

            if (existing >= 0) {
                id = existing;
            } else {
                if (id < 0) {
                    id = static_cast<ID>(_data->lines.size());
                    _data->lines.reserve(static_cast<std::size_t>(id + 128));
                    _data->lines.resize(static_cast<std::size_t>(id + 1));
                } else if (static_cast<std::size_t>(id) >= _data->lines.size()) {
                    _data->lines.resize(static_cast<std::size_t>(id + 1));
                }
                _data->lines[static_cast<std::size_t>(id)] = ln;
            }
            // Add to the currently active line set (activated in register_surface)
            _data->line_sets.add(id);
            return;
        }

        logging::error(false, "surface with id=", surface_id, " has not been defined for element with id=", element_id);
    }

    // allow negative id = automatically assign id
    if (id < 0) {
        id = _data->surfaces.size();
        _data->surfaces.reserve(id + 128);
        _data->surfaces.resize(id + 1);
    }
    logging::error(_data->surfaces[id] == nullptr, "surface with id=", id, " has already been defined");

    _data->surfaces[id] = surfptr;
    _data->surface_sets.add(id);
}

inline void Model::set_surface(const std::string& elset, ID surface_id) {
    logging::error(_data->elem_sets.has(elset), "element set with name=", elset, " has not been defined");

    for (const auto& el_id : *_data->elem_sets.get(elset)) {
        set_surface(-1, el_id, surface_id);
    }
}

template<typename T, typename... Args>
inline void Model::add_coordinate_system(const std::string& name, Args&&... args) {
    logging::error(!_data->coordinate_systems.has(name), "coordinate system with name=", name, " has already been defined");
    _data->coordinate_systems.activate<T>(name, args...);
}
