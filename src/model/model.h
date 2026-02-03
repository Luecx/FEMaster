#pragma once

#include "../bc/amplitude.h"
#include "../constraints/connector.h"
#include "../constraints/coupling.h"
#include "../cos/coordinate_system.h"
#include "geometry/surface/surface.h"
#include "./model_data.h"
#include "../constraints/constraint_groups.h"
#include "../core/types_cls.h"
#include "element/element.h"
#include "element/element_structural.h"

namespace fem {

namespace model{

struct Model {
    ModelDataPtr _data;

    // error out if not all elements have the same dimension (e.g. cannot use 1d, 2d and 3d elements at the same time)
    Dim element_dims = 0;

    // constructor which defines max elements and max nodes
    Model(ID max_nodes, ID max_elems, ID max_surfaces, ID max_integration_points = 0) :
        _data(std::make_shared<ModelData>(max_nodes, max_elems, max_surfaces, max_integration_points)){
        auto positions = _data->create_field("POSITION", FieldDomain::NODE, 6, false);
        positions->set_zero();
        _data->positions = positions;
    }

    // adding nodes and elements
    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);
    inline void set_surface(ID id, ID element_id, ID surface_id);
    inline void set_surface(const std::string& elset, ID surface_id);
    template<typename T, typename... Args>
    inline void set_beam_element(ID id, ID orientation_node, Args&&... args);

    // setting coordinate systems
    template<typename T, typename... Args>
    inline void add_coordinate_system(const std::string& name, Args&&... args);

    // add _couplings for structural problems
    void add_connector(const std::string& set1      , const std::string& set2, const std::string& coordinate_system, constraint::ConnectorType type);
    void add_coupling (const std::string& master_set, const std::string& slave_set, Dofs coupled_dofs, constraint::CouplingType type, bool is_surface);
    void add_tie      (const std::string& master_set, const std::string& slave_set, Precision distance, bool adjust);

    // -------------------- STRUCTURAL --------------------
    // load management
    void add_cload      (const std::string& nset, Vec6 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_cload      (ID id, Vec6 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_dload      (const std::string& sfset, Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_dload      (ID id, Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_pload      (const std::string& sfset, Precision load, const std::string& amplitude = "");
    void add_pload      (ID id, Precision load, const std::string& amplitude = "");
    void add_vload      (const std::string& elset, Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_vload      (ID id, Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_tload      (std::string& temp_field, Precision ref_temp);

    void define_amplitude(const std::string& name, bc::Interpolation interpolation);
    void add_amplitude_sample(const std::string& name, Precision time, Precision value);

    // support managment
    void add_support    (const std::string& nset, StaticVector<6> constraint, const std::string& orientation="");
    void add_support    (ID id, StaticVector<6> constraint, const std::string& orientation="");

    // filling in fields
    void set_field_temperature(const std::string& name, ID id, Precision value);

    // connecting materials with elements
    void solid_section(const std::string& set, const std::string& material);
    void beam_section (const std::string& set, const std::string& material, const std::string& profile, Vec3 orientation);
    void shell_section(const std::string& set, const std::string& material, Precision thickness);

    // features
    void add_point_mass_feature(const std::string& nset,
                                Precision mass,
                                Vec3 rotary_inertia,
                                Vec3 spring_constants,
                                Vec3 rotary_spring_constants);

    // stream output to console
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);

    // assigns sections to each element
    void assign_sections();

    // solving the given problem  set
    SystemDofIds  build_unconstrained_index_matrix();
    Field         build_integration_point_numeration();

    // building loads for every node including non existing ones
    Field                 build_load_matrix(std::vector<std::string> load_sets = {}, Precision time = 0);
    constraint::ConstraintGroups collect_constraints(SystemDofIds& system_dof_ids, const std::vector<std::string>& supp_sets = {});
    constraint::Equations build_constraints(SystemDofIds& system_dof_ids, std::vector<std::string> supp_sets = {});

    // matrices
    SparseMatrix  build_stiffness_matrix     (SystemDofIds& indices, const Field* stiffness_scalar = nullptr);
    SparseMatrix  build_geom_stiffness_matrix(SystemDofIds& indices,
                                              const Field& ip_stress,
                                              const Field* stiffness_scalar = nullptr);
    SparseMatrix  build_lumped_mass_matrix  (SystemDofIds& indices);


    std::tuple<Field, Field>       compute_ip_stress_strain(Field& displacement);
    std::tuple<Field, Field>       compute_stress_strain(Field& displacement);
    Field                          compute_compliance   (Field& displacement);
    Field                          compute_compliance_angle_derivative(Field& displacement);
    Field                          compute_volumes      ();
    DynamicMatrix                  compute_section_forces(Field& displacement);
};

#include "model.ipp"
} }    // namespace fem
