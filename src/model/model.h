#pragma once

#include "../bc/amplitude.h"
#include "../constraints/types/connector.h"
#include "../constraints/types/coupling.h"
#include "../cos/coordinate_system.h"
#include "geometry/surface/surface.h"
#include "./model_data.h"
#include "../constraints/constraint_groups.h"
#include "../core/types_cls.h"
#include "element/element.h"
#include "element/element_structural.h"

#include <utility>

namespace fem {
namespace model{
struct Model {
    ModelDataPtr _data;

    // error out if not all elements have the same dimension (e.g. cannot use 1d, 2d and 3d elements at the same time)
    Dim element_dims = 0;

    // constructor which defines max elements and max nodes
    Model(ID max_nodes, ID max_elems, ID max_surfaces, ID max_integration_points = 0) :
        _data(std::make_shared<ModelData>(max_nodes, max_elems, max_surfaces, max_integration_points)) {
        auto positions           = _data->create_field("POSITION"          , FieldDomain::NODE, 6, false);
        auto positions_reference = _data->create_field("POSITION_REFERENCE", FieldDomain::NODE, 6, false);
        positions->set_zero();
        positions_reference->set_zero();
        _data->positions           = positions;
        _data->positions_reference = positions_reference;
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
    void add_contact  (const std::string& master_set,
                       const std::string& slave_set,
                       Precision distance,
                       Precision penalty,
                       Precision clearance,
                       bool flip_normal);
    // rbm constraints: remove rigid-body motion on an element set
    void add_rbm      (const std::string& set);

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

    // inertia load (rigid-body field): target element set, optionally plus point-mass features
    void add_inertialload(const std::string& elset,
                          Vec3 center,
                          Vec3 center_acceleration,
                          Vec3 angular_velocity,
                          Vec3 angular_acceleration,
                          bool consider_point_masses = false);

    void define_amplitude(const std::string& name, bc::Interpolation interpolation);
    void add_amplitude_sample(const std::string& name, Precision time, Precision value);

    // support managment
    void add_support    (const std::string& nset, StaticVector<6> constraint, const std::string& orientation="");
    void add_support    (ID id, StaticVector<6> constraint, const std::string& orientation="");

    // connecting materials with elements
    void solid_section(const std::string& set, const std::string& material, const std::string& orientation = "");
    void beam_section (const std::string& set, const std::string& material, const std::string& profile, Vec3 orientation);
    void truss_section(const std::string& set, const std::string& material, Precision area);
    void shell_section(const std::string& set, const std::string& material, Precision thickness, const std::string& orientation = "");
    void shell_section_abd(const std::string& set,
                           const std::string& material,
                           Precision          thickness,
                           const Mat6&        abd,
                           const Mat2&        shear,
                           const std::string& orientation = "");

    // features
    void add_point_mass_feature(const std::string& nset,
                                Precision mass,
                                Vec3 rotary_inertia,
                                Vec3 spring_constants,
                                Vec3 rotary_spring_constants);

    // stream output to console
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);

    // Step-local element lifecycle.
    void step_begin();
    void step_end();

    // assigns sections to each element
    void assign_sections();

    // solving the given problem  set
    SystemDofIds  build_unconstrained_index_matrix();

    // building loads for every node including non existing ones
    Field                 build_load_matrix(std::vector<std::string> load_sets = {}, Precision time = 0);
    std::vector<std::pair<bc::Amplitude::Ptr, Field>> build_load_basis(std::vector<std::string> load_sets = {});
    constraint::ConstraintGroups collect_constraints(SystemDofIds& system_dof_ids, const std::vector<std::string>& supp_sets = {});
    constraint::Equations build_constraints(SystemDofIds& system_dof_ids, std::vector<std::string> supp_sets = {});

    // matrices
    SparseMatrix  build_stiffness_matrix     (SystemDofIds& indices, const Field* stiffness_scalar = nullptr);
    SparseMatrix  build_tangent_stiffness_matrix(SystemDofIds& indices,
                                                  NodeData& nodal_forces,
                                                  const Field& displacement,
                                                  const Field* stiffness_scalar = nullptr);
    SparseMatrix  build_geom_stiffness_matrix(SystemDofIds& indices,
                                              const Field& ip_stress,
                                              const Field* stiffness_scalar = nullptr);
    SparseMatrix  build_lumped_mass_matrix  (SystemDofIds& indices);
    Field         build_internal_force_nonlinear(const Field& ip_stress);

    std::tuple<Field, Field>       compute_stress_ip(Field& displacement, bool use_green_lagrange_nl = false);
    Field                          compute_stress_state(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field>       compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field>       compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field>       compute_ip_stress_strain(Field& displacement);
    Field                          compute_ip_stress_nonlinear(Field& displacement);
    std::tuple<Field, Field>       compute_stress_strain(Field& displacement);
    std::tuple<Field, Field>       compute_shell_stress_surfaces(Field& displacement);
    Field                          compute_shell_resultants(Field& displacement);
    Field                          compute_compliance   (Field& displacement);
    Field                          compute_compliance_angle_derivative(Field& displacement);
    Field                          compute_volumes      ();
    Field                          compute_section_forces(Field& displacement);
    Field                          compute_shear_flow(Field& displacement);
};

#include "model.ipp"
} }    // namespace fem
