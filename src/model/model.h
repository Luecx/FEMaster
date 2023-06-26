#pragma once

#include "../core/core.h"
#include "element.h"
#include "sets.h"

namespace fem {

namespace model{

struct Model {
    const ID max_nodes;
    const ID max_elements;

    // nodes and elements
    NodeData node_coords;
    std::vector<ElementInterface*> elements{};

    // sets to group nodes and elements
    Sets node_sets{"NALL"};
    Sets elem_sets{"EALL"};

    // error out if not all elements have the same dimension (e.g. cannot use 1d, 2d and 3d elements at the same time)
    Dim element_dims = 0;

    // constructor which defines max elements and max nodes
    Model(ID max_nodes, ID max_elems) : max_nodes(max_nodes), max_elements(max_elems), node_coords(max_nodes, 3){
        elements.resize(max_elems);
    }

    // adding nodes and elements
    void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    void set_element(ID id, Args&&... args);


    // functions to generate the stiffness matrix as well as the boundary conditions

    SparseMatrix stiffness(DynamicVector &constraints){
        SparseMatrixBuilder tripplets{};
    }
    DynamicVector load(DynamicVector &constraints, DynamicVector &load){
        return DynamicVector{1};
    }

};

#include "model.ipp"

}

}    // namespace fem

