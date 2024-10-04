inline void Model::set_node(ID id, Precision x, Precision y, Precision z) {
    logging::error(id < max_nodes, "internal error; allocated less data than required. id=", id, " exceeds maximum limit");
    logging::error(node_coords(id,0) == 0 &&
                  node_coords(id,1) == 0 &&
                  node_coords(id,2) == 0, "node with id=", id, " seems to define non-zero values twice");
    node_coords(id, 0) = x;
    node_coords(id, 1) = y;
    node_coords(id, 2) = z;
    node_sets.add(id);
}

template<typename T, typename... Args>
inline void Model::set_element(ID id, Args&&... args) {
    logging::error(id < max_elements,
              "internal error; allocated less data than required. id=", id, " exceeds maximum limit");
    auto el = ElementPtr {new T {id, {args...}}};
    if (element_dims == 0) {
        element_dims = el->dimensions();
    }
    logging::error(element_dims == el->dimensions(), "different number of dimensions across elements of model; element with id=",
              id, " has ", el->dimensions(), " while all other elements have ", element_dims, " dimensions");
    logging::error(this->elements[id] == nullptr, "element with id=", id, " has already been defined");

    this->elements[id] = el;
    this->elem_sets.add(id);
}

inline void Model::set_surface(ID id, ID element_id, ID surface_id) {
    logging::error(id < max_surfaces,
                   "internal error; allocated less data than required. id=", id, " exceeds maximum limit");

    auto elptr = elements[element_id];
    logging::error(elptr != nullptr, "element with id=", element_id, " has not been defined");
    auto surfptr = elptr->surface(surface_id);
    logging::error(surfptr != nullptr, "surface with id=", surface_id, " has not been defined for element with id=", element_id);
    logging::error(this->surfaces[id] == nullptr, "surface with id=", id, " has already been defined");

    this->surfaces[id] = surfptr;
    this->surface_sets.add(id);
}