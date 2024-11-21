#pragma once

#include "../constraints/connector.h"
#include "../constraints/coupling.h"
#include "../constraints/tie.h"
#include "../core/core.h"
#include "../cos/coordinate_system.h"
#include "element/element.h"
#include "element/element_structural.h"
#include "geometry/surface/surface.h"
#include "../data/dict.h"
#include "../data/region.h"
#include "../data/sets.h"

#include "../data/elem_data_dict.h"
#include "../data/node_data_dict.h"
#include "./model_data.h"

namespace fem {

namespace model{

struct Model {
    ModelDataPtr _data;

    // constraints
    std::vector<constraint::Connector>       _connectors {};
    std::vector<constraint::Coupling>        _couplings {};
    std::vector<constraint::Tie>             _ties {};

    // error out if not all elements have the same dimension (e.g. cannot use 1d, 2d and 3d elements at the same time)
    Dim element_dims = 0;

    // storage for all the sets
    Dict<NodeData>  _load_sets   {};
    Dict<NodeData>  _support_sets{};

    // storage for other fields
    Dict<NodeData>   _fields_temperature;

    // constructor which defines max elements and max nodes
    Model(ID max_nodes, ID max_elems, ID max_surfaces) :
        _data(std::make_shared<ModelData>(max_nodes, max_elems, max_surfaces)){

        // initialize the node data
        _data->node_data.create(NodeDataEntries::POSITION, 6);
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

    // add _couplings for structural problems
    void add_connector(const std::string& set1, const std::string& set2, const std::string& coordinate_system, constraint::ConnectorType type);
    void add_coupling(const std::string& master_set, const std::string& slave_set, Dofs coupled_dofs, constraint::CouplingType type);
    void add_tie(const std::string& master_set, const std::string& slave_set, Precision distance, bool adjust);

    // -------------------- STRUCTURAL --------------------
    // load management
    void add_cload      (const std::string& nset, Vec6 load);
    void add_cload      (ID id, Vec6 load);
    void add_dload      (const std::string& sfset, Vec3 load);
    void add_dload      (ID id, Vec3 load);
    void add_vload      (const std::string& elset, Vec3 load);
    void add_vload      (ID id, Vec3 load);
    void add_tload      (std::string& temp_field, Precision ref_temp);

    // support managment
    void add_support    (const std::string& nset, StaticVector<6> constraint);
    void add_support    (ID id, StaticVector<6> constraint);

    // filling in fields
    void set_field_temperature(const std::string& name, ID id, Precision value);

    // connecting materials with elements
    void solid_section(const std::string& set, const std::string& material);

    // stream output to console
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);

    // assigns sections to each element
    void assign_sections();

    // solving the given problem  set
    SystemDofIds  build_unconstrained_index_matrix();

    // building constraints and loads for every node including non existing ones
    NodeData    build_support_matrix (std::vector<std::string> support_sets = {});
    NodeData    build_load_matrix    (std::vector<std::string> load_sets = {});

    // matrices
    SparseMatrix  build_constraint_matrix   (SystemDofIds& indices, Precision characteristic_stiffness=1.0);
    SparseMatrix  build_stiffness_matrix    (SystemDofIds& indices, ElementData stiffness_scalar = ElementData(0,0));
    SparseMatrix  build_lumped_mass_matrix  (SystemDofIds& indices);

    std::tuple<NodeData, NodeData> compute_stress_strain(NodeData& displacement);
    ElementData                    compute_compliance   (NodeData& displacement);
    ElementData                    compute_compliance_angle_derivative(NodeData& displacement);
    ElementData                    compute_volumes      ();
};

#include "model.ipp"
} }    // namespace fem

