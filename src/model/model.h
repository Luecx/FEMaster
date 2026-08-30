/**
 * @file model.h
 * @brief Declares the high-level FEM model interface.
 */

#pragma once

#include "../bc/amplitude.h"
#include "../bc/dirichlet/dirichlet.h"
#include "../bc/load.h"
#include "../bc/robin/robin.h"
#include "../constraints/constraint_groups.h"
#include "../core/types_cls.h"
#include "../cos/coordinate_system.h"
#include "element/element.h"
#include "element/element_structural.h"
#include "geometry/surface/surface.h"
#include "instance.h"
#include "model_data.h"
#include "part.h"

#include <ostream>
#include <utility>

namespace fem::model {

/** @brief Nodal thermal RHS and symbolic Robin equation contributions. */
struct ThermalLoadData {
    Field               rhs;
    bc::RobinEquations  equations;
};

struct Model {
    static constexpr const char* DEFAULT_PART_NAME     = "__DEFAULT_PART__";
    static constexpr const char* DEFAULT_INSTANCE_NAME = "__DEFAULT_INSTANCE__";

    ModelDataPtr _data;

    Model();
    Model(const Model&)            = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) noexcept        = default;

    void add_part(const std::string& name);
    void add_instance(const std::string& name,
                      const std::string& part,
                      Vec3 translation = Vec3::Zero(),
                      Mat3 rotation = Mat3::Identity());
    void compile();

    ID compiled_node_id   (ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_element_id(ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_surface_id(ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_line_id   (ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;

    ID compiled_node_id   (const std::string& reference) const;
    ID compiled_element_id(const std::string& reference) const;
    ID compiled_surface_id(const std::string& reference) const;

    void add_nodes_to_part_set        (const std::string& set, const std::string& source);
    void add_nodes_to_assembly_set    (const std::string& set, const std::string& source);
    void add_elements_to_part_set     (const std::string& set, const std::string& source);
    void add_elements_to_assembly_set (const std::string& set, const std::string& source);
    void add_surfaces_to_part_set     (const std::string& set, const std::string& source);
    void add_surfaces_to_assembly_set (const std::string& set, const std::string& source);

    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);
    inline void set_surface(ID id, ID element_id, ID surface_id);
    inline void set_surface(const std::string& elset, ID surface_id);

    void add_coordinate_system(cos::CoordinateSystem::Ptr coordinate_system);
    void add_material(material::Material::Ptr material);
    void add_profile(Profile::Ptr profile);
    void add_section(Section::Ptr section);

    void add_load     (bc::Load::Ptr load);
    void add_amplitude(bc::Amplitude::Ptr amplitude);
    void add_support  (bc::Dirichlet::Ptr support);

    void step_begin();
    void step_end();
    void assign_sections();
    void build_shell_element_normals(Precision equalize_angle_degrees = Precision(20));

    [[nodiscard]] Index maximum_material_state_size() const;
    void initialize_material_state(Field& state) const;

    SystemDofIds build_unconstrained_index_matrix();
    SystemDofIds build_thermal_index_matrix();
    Field build_field_mapping_weights(bool structural = true, bool thermal = true) const;

    Field build_structural_load_matrix(
        const std::vector<std::string>& load_sets = {},
        Precision time = 0);
    ThermalLoadData build_thermal_loads(
        const std::vector<std::string>& load_sets,
        Precision time = 0);

    constraint::ConstraintGroups collect_constraints(
        SystemDofIds& system_dof_ids,
        const std::vector<std::string>& supp_sets = {});
    constraint::ConstraintGroups collect_temperature_constraints(
        const std::vector<std::string>& supp_sets = {});

    std::vector<std::pair<bc::Amplitude::Ptr, Field>> build_load_basis(
        std::vector<std::string> load_sets = {});

    SparseMatrix build_stiffness_matrix(
        SystemDofIds& indices,
        const Field* stiffness_scalar = nullptr);
    SparseMatrix build_conductivity_matrix(SystemDofIds& indices);
    SparseMatrix build_thermal_matrix(
        SystemDofIds& indices,
        const bc::RobinEquations& equations);
    SparseMatrix build_tangent_stiffness_matrix(
        SystemDofIds& indices,
        NodeData& nodal_forces,
        const Field& displacement,
        const Field* stiffness_scalar = nullptr);
    SparseMatrix build_geom_stiffness_matrix(
        SystemDofIds& indices,
        const Field& ip_stress,
        const Field* stiffness_scalar = nullptr);
    void build_internal_force_nonlinear(
        SystemDofIds& indices,
        NodeData& nodal_forces,
        const Field& displacement);
    SparseMatrix build_lumped_mass_matrix(SystemDofIds& indices);

    Field compute_stress_state(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl = false);
    Field compute_shell_resultants(Field& displacement);
    Field compute_compliance(Field& displacement);
    Field compute_compliance_angle_derivative(Field& displacement);
    Field compute_volumes();
    Field compute_section_forces(Field& displacement);
    Field compute_shear_flow(Field& displacement);
    Field compute_heat_flux(const Field& temperature);

    void print_overview() const;
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);
};

#include "model.ipp"

} // namespace fem::model
