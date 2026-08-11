/**
 * @file model.h
 * @brief Declares the high-level finite-element model facade.
 *
 * `Model` builds and operates on the shared topology, fields, materials,
 * sections, loads and constraints stored by `ModelData`. It provides model
 * construction, element lifecycle and global assembly entry points used by load
 * cases, while individual elements, sections and constraints retain their local
 * mechanical responsibilities.
 *
 * Material-point enumeration belongs to `ModelData`. This facade determines the
 * state width required by assigned constitutive laws and initializes every
 * enumerated material-point row; nonlinear trial ownership remains with the
 * load-case state manager. Surface-to-surface contact definitions are created
 * here and assembled by `constraint::Contact` during nonlinear analysis.
 *
 * @see ModelData
 * @see Element
 * @see constraint::Contact
 * @see loadcase::tools::NonlinearStateManager
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#pragma once

#include "../bc/amplitude.h"
#include "../constraints/constraint_groups.h"
#include "../constraints/types/connector.h"
#include "../constraints/types/coupling.h"
#include "../core/types_cls.h"
#include "../cos/coordinate_system.h"
#include "element/element.h"
#include "element/element_structural.h"
#include "geometry/surface/surface.h"
#include "model_data.h"

#include <utility>

namespace fem {
namespace model {

/**
 * @brief High-level owner and assembly interface for one finite-element model.
 *
 * Persistent entities and fields are held in the shared `ModelData` repository
 * so elements, load cases, constraints and writers observe the same model state.
 * The facade resolves named regions and model definitions, coordinates element
 * step lifecycle, and exposes global load, constraint, matrix and result assembly.
 *
 * Element-local kinematics and constitutive operations remain responsibilities of
 * the corresponding element/section implementations. Material and contact trial
 * histories are coordinated by the nonlinear load-case state manager rather than
 * being owned directly by this facade.
 */
struct Model {
    // Shared persistent model repository and common element dimensionality
    ModelDataPtr _data;
    Dim          element_dims = 0;

    // Construction with fixed topology capacities and initialized position fields
    Model(ID max_nodes, ID max_elems, ID max_surfaces, ID max_integration_points = 0) :
        _data(std::make_shared<ModelData>(max_nodes, max_elems, max_surfaces, max_integration_points)) {
        auto positions           = _data->create_field("POSITION"          , FieldDomain::NODE, 6, false);
        auto positions_reference = _data->create_field("POSITION_REFERENCE", FieldDomain::NODE, 6, false);
        positions->set_zero();
        positions_reference->set_zero();
        _data->positions           = positions;
        _data->positions_reference = positions_reference;
    }

    Model(const Model&)            = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) noexcept        = default;

    // Topology construction
    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);
    inline void set_surface(ID id, ID element_id, ID surface_id);
    inline void set_surface(const std::string& elset, ID surface_id);
    template<typename T, typename... Args>
    inline void set_beam_element(ID id, ID orientation_node, Args&&... args);

    // Coordinate-system construction
    template<typename T, typename... Args>
    inline void add_coordinate_system(const std::string& name, Args&&... args);

    // Structural connectors, couplings, ties and unilateral surface contact
    void add_connector(
        const std::string& set1,
        const std::string& set2,
        const std::string& coordinate_system,
        constraint::ConnectorType type);
    void add_coupling(
        const std::string& master_set,
        const std::string& slave_set,
        Dofs coupled_dofs,
        constraint::CouplingType type,
        bool is_surface);
    void add_tie(
        const std::string& master_set,
        const std::string& slave_set,
        Precision distance,
        bool adjust);
    void add_contact(
        const std::string& master_set,
        const std::string& slave_set,
        Precision penalty,
        Precision clearance,
        bool flip_normal);

    // Rigid-body-motion suppression on an element set
    void add_rbm(const std::string& set);

    // Structural nodal, surface, volume and temperature loads
    void add_cload(const std::string& nset , Vec6 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_cload(ID id                   , Vec6 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_dload(const std::string& sfset, Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_dload(ID id                   , Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_pload(const std::string& sfset, Precision load, const std::string& amplitude = "");
    void add_pload(ID id                   , Precision load, const std::string& amplitude = "");
    void add_vload(const std::string& elset, Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_vload(ID id                   , Vec3 load, const std::string& orientation = "", const std::string& amplitude = "");
    void add_tload(std::string& temp_field , Precision ref_temp);
    void add_inertialload(const std::string& elset,
                          Vec3 center,
                          Vec3 center_acceleration,
                          Vec3 angular_velocity,
                          Vec3 angular_acceleration,
                          bool consider_point_masses = false);

    // Time-dependent load amplitudes
    void define_amplitude(const std::string& name, bc::Interpolation interpolation);
    void add_amplitude_sample(const std::string& name, Precision time, Precision value);

    // Nodal supports
    void add_support(const std::string& nset, StaticVector<6> constraint, const std::string& orientation = "");
    void add_support(ID id, StaticVector<6> constraint, const std::string& orientation = "");

    // Element section assignments
    void solid_section    (const std::string& set, const std::string& material, const std::string& orientation = "");
    void beam_section     (const std::string& set, const std::string& material, const std::string& profile, Vec3 orientation);
    void truss_section    (const std::string& set, const std::string& material, Precision area);
    void shell_section    (const std::string& set, const std::string& material, Precision thickness, const std::string& orientation = "", Index csys_axis = 0);
    void shell_section_abd(const std::string& set, const std::string& material, Precision thickness, const Mat6& abd,
                           const Mat2& shear, const std::string& orientation = "", Index csys_axis = 0);

    // Point-based mass and spring features
    void add_point_mass_feature(const std::string& nset,
                                Precision mass,
                                Vec3 rotary_inertia,
                                Vec3 spring_constants,
                                Vec3 rotary_spring_constants);

    // Human-readable model summary
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);

    // Step-local element lifecycle and section assignment
    void step_begin();
    void step_end();
    void assign_sections();

    // Smoothed element-nodal reference normals for shell directors
    void build_shell_element_normals(Precision equalize_angle_degrees = Precision(20));

    // Material-point state sizing and reference-history initialization. Rows use
    // the global element/IP/material-point enumeration stored by ModelData.
    [[nodiscard]] Index maximum_material_state_size() const;
    void initialize_material_state(Field& state) const;

    // Active global DOFs, loads and linear constraint equations
    SystemDofIds build_unconstrained_index_matrix();
    Field build_load_matrix(
        std::vector<std::string> load_sets = {}, Precision time = 0);
    std::vector<std::pair<bc::Amplitude::Ptr, Field>> build_load_basis(
        std::vector<std::string> load_sets = {});
    constraint::ConstraintGroups collect_constraints(
        SystemDofIds& system_dof_ids,
        const std::vector<std::string>& supp_sets = {});
    constraint::Equations build_constraints(SystemDofIds& system_dof_ids,
                                             std::vector<std::string> supp_sets = {});

    // Global structural/contact matrix and nonlinear internal-force assembly
    SparseMatrix build_stiffness_matrix(
        SystemDofIds& indices,
        const Field* stiffness_scalar = nullptr);
    SparseMatrix build_tangent_stiffness_matrix(
        SystemDofIds& indices,
        NodeData& nodal_forces,
        const Field& displacement,
        const Field* stiffness_scalar = nullptr,
        bool assemble_contact = true);
    SparseMatrix build_geom_stiffness_matrix(
        SystemDofIds& indices,
        const Field& ip_stress,
        const Field* stiffness_scalar = nullptr);
    void build_internal_force_nonlinear(
        SystemDofIds& indices,
        NodeData& nodal_forces,
        const Field& displacement,
        bool assemble_contact = true);
    SparseMatrix build_lumped_mass_matrix(
        SystemDofIds& indices);

    // Structural result recovery and derived element/nodal fields
    Field                    compute_stress_state(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl = false);
    Field                    compute_shell_resultants(Field& displacement);
    Field                    compute_compliance(Field& displacement);
    Field                    compute_compliance_angle_derivative(Field& displacement);
    Field                    compute_volumes();
    Field                    compute_section_forces(Field& displacement);
    Field                    compute_shear_flow(Field& displacement);
};

#include "model.ipp"

} // namespace model
} // namespace fem
