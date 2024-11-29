//
// Created by f_eggers on 11.11.2024.
//

#ifndef MODEL_DATA_H
#define MODEL_DATA_H

#include "../bc/load.h"
#include "../bc/load_collector.h"
#include "../bc/support.h"
#include "../bc/support_collector.h"
#include "../constraints/connector.h"
#include "../constraints/coupling.h"
#include "../constraints/equation.h"
#include "../constraints/tie.h"
#include "../core/types.h"
#include "../data/dict.h"
#include "../data/elem_data_dict.h"
#include "../data/node_data_dict.h"
#include "../data/region.h"
#include "../data/sets.h"
#include "../section/section.h"

namespace fem::model {

struct ModelData {

    // constants
    ID max_nodes;
    ID max_elems;
    ID max_surfaces;

    // geometric entities
    std::vector<ElementPtr> elements;
    std::vector<SurfacePtr> surfaces;

    // sections with their respective data
    std::vector<Section::Ptr> sections;

    // data for nodes, elements. Only one allowed per field
    // other data is stored in sets
    NodeDataDict node_data;
    ElemDataDict elem_data;

    // field data
    NodeFieldDict node_fields;
    ElemFieldDict elem_fields;

    // regions
    Sets<NodeRegion>     node_sets    {SET_NODE_ALL};
    Sets<ElementRegion>  elem_sets    {SET_ELEM_ALL};
    Sets<SurfaceRegion>  surface_sets {SET_SURF_ALL};

    // other things which are named
    Dict<material::Material> materials;
    Dict<cos::CoordinateSystem> coordinate_systems;

    // constraints
    std::vector<constraint::Connector>  connectors {};
    std::vector<constraint::Coupling>   couplings {};
    std::vector<constraint::Tie>        ties {};
    std::vector<constraint::Equation>   equations {};

    // mechanical loads and supports
    Sets<SupportCollector> supp_cols{};
    Sets<LoadCollector>    load_cols{};
    // Dict<NodeData>  load_cols{};

    //
    // Constructor
    //
    ModelData(ID max_nodes, ID max_elems, ID max_surfaces)
        : max_nodes(max_nodes),
          max_elems(max_elems),
          max_surfaces(max_surfaces),
            node_data(max_nodes),
            elem_data(max_elems) {
        elements.resize(max_elems);
        surfaces.resize(max_surfaces);
    }

    // managing of data for nodes and elements (not for fields)
    void create_data(const NodeDataEntries key, int entries) {
        node_data.create(key, entries);
    }
    void create_data(const ElementDataEntries key, int entries) {
        elem_data.create(key, entries);
    }

    NodeData& get(const NodeDataEntries key) {
        return node_data.get(key);
    }
    ElementData& get(const ElementDataEntries key) {
        return elem_data.get(key);
    }

    void remove(const NodeDataEntries key) {
        node_data.remove(key);
    }
    void remove(const ElementDataEntries key) {
        elem_data.remove(key);
    }

};

}


#endif //MODEL_DATA_H
