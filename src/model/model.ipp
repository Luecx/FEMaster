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

inline void Model::set_surface(ID id, ID element_id, ID surface_id) {
    logging::error(id < _data->max_surfaces, "internal error; allocated less data than required. id=", id, " exceeds maximum limit");

    auto& elptr = _data->elements[element_id];
    logging::error(elptr != nullptr, "element with id=", element_id, " has not been defined");
    auto surfptr = elptr->surface(surface_id);
    logging::error(surfptr != nullptr, "surface with id=", surface_id, " has not been defined for element with id=", element_id);

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
