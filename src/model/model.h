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
    Sets<std::vector<ID>> node_sets{"NALL"};
    Sets<std::vector<ID>> elem_sets{"EALL"};

    // error out if not all elements have the same dimension (e.g. cannot use 1d, 2d and 3d elements at the same time)
    Dim element_dims = 0;

    // storage for all the load collectors
    Sets<NodeData> load_collectors{"LALL"};
    Sets<NodeData> support_collectors{"SALL"};

    // constructor which defines max elements and max nodes
    Model(ID max_nodes, ID max_elems) : max_nodes(max_nodes), max_elements(max_elems), node_coords(max_nodes, 3){
        elements.resize(max_elems);
    }

    // adding nodes and elements
    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);
    inline void activate_node_set(const std::string &name);
    inline void activate_element_set(const std::string &name);
    inline void activate_load_set(const std::string &name);
    inline void activate_support_set(const std::string &name);

    // solving the given problem set
    IndexMatrix build_dof_index_vector(DynamicVector& constraints);
    DynamicVector solve(DynamicVector &constraints, DynamicVector &loads);
    SparseMatrix stiffness();
    void apply_constraints(DynamicVector& constraints, SparseMatrix& matrix, DynamicVector& loads);
};



#include "model.ipp"

}

}    // namespace fem

