#pragma once

#include "../core/core.h"
#include "element.h"
#include "sets.h"
#include "../constraints/coupling.h"

namespace fem {

namespace model{

struct Model {
    const ID max_nodes;
    const ID max_elements;

    // nodes and elements
    NodeData node_coords;
    std::vector<ElementPtr> elements{};
    std::vector<Coupling>   couplings{};

    // sets to group nodes and elements
    Sets<std::vector<ID>> node_sets{"NALL"};
    Sets<std::vector<ID>> elem_sets{"EALL"};

    // error out if not all elements have the same dimension (e.g. cannot use 1d, 2d and 3d elements at the same time)
    Dim element_dims = 0;

    // storage for all the load collectors
    Sets<NodeData> load_sets;
    Sets<NodeData> support_sets;
    Sets<material::Material> materials;

    // constructor which defines max elements and max nodes
    Model(ID max_nodes, ID max_elems) :
        max_nodes   (max_nodes),
        max_elements(max_elems),
        node_coords (max_nodes, 3),
        load_sets   ("LALL", max_nodes, 6),
        support_sets("SALL", max_nodes, 6){

        // clear node coords
        node_coords.setZero();

        // clear sets to make sure
        load_sets   .current().setZero();
        support_sets.current().setZero();
        support_sets.current().fill(std::numeric_limits<Precision>::quiet_NaN());

        // resize elements to avoid later reallocation
        elements.resize(max_elems);
    }

    // adding nodes and elements
    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);

    // add couplings
    void add_coupling(const std::string& master_set, const std::string& slave_set, Dofs coupled_dofs, CouplingType type);

    // set management
    void activate_node_set   (const std::string &name);
    void activate_element_set(const std::string &name);
    void activate_load_set   (const std::string &name);
    void activate_support_set(const std::string &name);
    void activate_material   (const std::string& name);

    // load management
    void add_cload      (const std::string& nset, StaticVector<3> load);
    void add_cload      (ID id, StaticVector<3> load);
    void add_vload      (const std::string& elset, StaticVector<3> load);
    void add_vload      (ID id, StaticVector<3> load);

    // support managment
    void add_support    (const std::string& nset, StaticVector<6> constraint);
    void add_support    (const std::string& nset, StaticVector<3> displacement);
    void add_support_rot(const std::string& nset, StaticVector<3> rotation);
    void add_support    (ID id, StaticVector<6> constraint);
    void add_support    (ID id, StaticVector<3> displacement);
    void add_support_rot(ID id, StaticVector<3> rotation);
    void add_support    (ID id, Dim dim, Precision displacement = 0);

    // access to active sets
    material::Material& active_material();
    NodeData&           active_loads();
    NodeData&           active_supports();
    std::vector<ID>&    active_nodeset();
    std::vector<ID>&    active_elemset();

    Sets<std::vector<ID>>&  nodesets();
    Sets<std::vector<ID>>&  elemsets();

    // connecting materials with elements
    void solid_section(const std::string& set, const std::string& material);

    // stream output to console
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);

    // solving the given problem  set
    SystemDofIds  build_unconstrained_index_matrix();

    // building constraints and loads for every node including non existing ones
    NodeData    build_support_matrix (std::vector<std::string> support_sets = {"SALL"});
    NodeData    build_load_matrix    (std::vector<std::string> load_sets = {"LALL"});

    // matrices
    SparseMatrix  build_constraint_matrix   (SystemDofIds& indices, Precision characteristic_stiffness=1.0);
    SparseMatrix  build_stiffness_matrix    (SystemDofIds& indices, ElementData stiffness_scalar = ElementData(0,0));
    SparseMatrix  build_lumped_mass_matrix  (SystemDofIds& indices);

    std::tuple<NodeData, NodeData> compute_stress_strain(NodeData& displacement);
    ElementData                    compute_compliance   (NodeData& displacement);
    ElementData                    compute_volumes      ();
};



#include "model.ipp"

}

}    // namespace fem

