#pragma once

#include "../core/core.h"
#include "element.h"
#include "surface/surface.h"
#include "sets.h"
#include "../constraints/coupling.h"
#include "../constraints/tie.h"
#include "../constraints/connector.h"
#include "../cos/coordinate_system.h"

namespace fem {

namespace model{

struct Model {
    const ID max_nodes;
    const ID max_elements;
    const ID max_surfaces;

    // nodes and elements
    NodeData node_coords;
    std::vector<ElementPtr> elements{};
    std::vector<SurfacePtr> surfaces{};


    // constraints
    std::vector<Connector>       connectors{};
    std::vector<Coupling>        couplings{};
    std::vector<constraint::Tie> ties{};

    // sets to group nodes and elements and everything that has a name
    Sets<std::vector<ID>> node_sets{SET_NODE_ALL};
    Sets<std::vector<ID>> elem_sets{SET_ELEM_ALL};
    Sets<std::vector<ID>> surface_sets{SET_SURF_ALL};

    Sets<cos::CoordinateSystemPtr> coordinate_systems;

    // error out if not all elements have the same dimension (e.g. cannot use 1d, 2d and 3d elements at the same time)
    Dim element_dims = 0;

    // storage for all the load collectors
    Sets<NodeData> load_sets;
    Sets<NodeData> support_sets;
    Sets<material::Material> materials;

    // constructor which defines max elements and max nodes
    Model(ID max_nodes, ID max_elems, ID max_surfaces) :
        max_nodes   (max_nodes),
        max_elements(max_elems),
        max_surfaces(max_surfaces),
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
        surfaces.resize(max_surfaces);
    }

    // adding nodes and elements
    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);
    inline void set_surface(ID id, ID element_id, ID surface_id);
    inline void set_surface(const std::string& elset, ID surface_id);

    // setting coordinate systems
    template<typename T, typename... Args>
    inline void add_coordinate_system(const std::string& name, Args&&... args);

    // add couplings
    void add_connector(const std::string& set1, const std::string& set2, const std::string& coordinate_system, ConnectorType type);
    void add_coupling(const std::string& master_set, const std::string& slave_set, Dofs coupled_dofs, CouplingType type);
    void add_tie(const std::string& master_set, const std::string& slave_set, Precision distance, bool adjust);

    // set management
    void activate_node_set   (const std::string &name);
    void activate_element_set(const std::string &name);
    void activate_surface_set(const std::string &name);
    void activate_load_set   (const std::string &name);
    void activate_support_set(const std::string &name);
    void activate_material   (const std::string& name);

    // load management
    void add_cload      (const std::string& nset, Vec3 load);
    void add_cload      (ID id, Vec3 load);
    void add_dload      (const std::string& sfset, Vec3 load);
    void add_dload      (ID id, Vec3 load);
    void add_vload      (const std::string& elset, Vec3 load);
    void add_vload      (ID id, Vec3 load);

    // support managment
    void add_support    (const std::string& nset, StaticVector<6> constraint);
    void add_support    (const std::string& nset, Vec3 displacement);
    void add_support_rot(const std::string& nset, Vec3 rotation);
    void add_support    (ID id, StaticVector<6> constraint);
    void add_support    (ID id, Vec3 displacement);
    void add_support_rot(ID id, Vec3 rotation);
    void add_support    (ID id, Dim dim, Precision displacement = 0);

    // access to active sets
    material::Material& active_material();
    NodeData&           active_loads();
    NodeData&           active_supports();
    std::vector<ID>&    active_nodeset();
    std::vector<ID>&    active_elemset();
    std::vector<ID>&    active_surfset();

    Sets<std::vector<ID>>&  nodesets();
    Sets<std::vector<ID>>&  elemsets();
    Sets<std::vector<ID>>&  surfsets();

    // connecting materials with elements
    void solid_section(const std::string& set, const std::string& material);

    // stream output to console
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);

    // solving the given problem  set
    SystemDofIds  build_unconstrained_index_matrix();

    // building constraints and loads for every node including non existing ones
    NodeData    build_support_matrix (std::vector<std::string> support_sets = {SET_SUPP_ALL});
    NodeData    build_load_matrix    (std::vector<std::string> load_sets = {SET_LOAD_ALL});

    // matrices
    SparseMatrix  build_constraint_matrix   (SystemDofIds& indices, Precision characteristic_stiffness=1.0);
    SparseMatrix  build_stiffness_matrix    (SystemDofIds& indices, ElementData stiffness_scalar = ElementData(0,0));
    SparseMatrix  build_lumped_mass_matrix  (SystemDofIds& indices);

    std::tuple<NodeData, NodeData> compute_stress_strain(NodeData& displacement);
    ElementData                    compute_compliance   (NodeData& displacement);
    ElementData                    compute_volumes      ();
};

#include "model.ipp"
} }    // namespace fem

