/**
 * @file model.h
 * @brief Declares the high-level FEM model interface.
 *
 * `ModelData` is the single owner of semantic definitions, part/instance
 * topology and dense solver storage. `Model` provides behavior around that data:
 * construction, the one-way part/instance compile pass, assembly and result
 * recovery. No parallel part, instance or compiled-state mirror is kept here.
 *
 * The `compiled` flag has one deliberately narrow meaning: the semantic
 * `Part`/`Instance` graph has already been flattened into the dense assembly and
 * may no longer change. Shared definition objects such as materials, profiles,
 * coordinate systems and compiled section assignments are not inherently frozen
 * by that transition and may be registered afterwards when their own
 * dependencies permit it.
 *
 * Constraint objects are intentionally not constructed or dispatched by Model.
 * Parsers and direct C++ callers create the concrete constraint type themselves
 * and append it to the corresponding collection in `ModelData`. This keeps the
 * model interface focused on model behavior rather than input-format factories.
 *
 * @see ModelData
 * @see Part
 * @see Instance
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../bc/amplitude.h"
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

/**
 * @brief Behavioral interface around semantic and compiled FEM model data.
 *
 * Before `compile()`, topology construction operates on the active Part stored
 * in `_data->parts`. The compile pass creates dense assembly nodes, elements,
 * surfaces and lines for every Instance and fills each instance-local-to-global
 * id map. Recompilation is intentionally unsupported because any later
 * assembly-level set, field, load or constraint may already reference those
 * dense identifiers.
 */
struct Model {
    // Root-level unqualified topology is stored as an implicit Part and embedded
    // through an identity Instance. This lets all semantic entity references use
    // one consistent namespace rule without a separate Assembly object.
    static constexpr const char* DEFAULT_PART_NAME     = "__DEFAULT_PART__";
    static constexpr const char* DEFAULT_INSTANCE_NAME = "__DEFAULT_INSTANCE__";

    // Single source of truth for both semantic model definitions and compiled
    // solver-facing storage.
    ModelDataPtr _data;

    // constructors
    Model();
    Model(const Model&)            = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) noexcept        = default;

    // adding parts and instances which can hold elements, nodes and sets.
    // compile() will later turn this into dense data inside _data.
    void add_part(const std::string& name);
    void add_instance(const std::string& name,
                      const std::string& part,
                      Vec3 translation = Vec3::Zero(),
                      Mat3 rotation = Mat3::Identity());


    // Flattens all Parts through their Instances into dense assembly data.
    // The operation may be called exactly once. It creates deterministic dense
    // identifiers, transforms nodal geometry and section orientations, rewires
    // copied element/surface/line connectivity, materializes inherited sets and
    // initializes element-nodal/IP/MP enumeration.
    void compile();

    // Resolve semantic local identifiers through the compile maps. Unqualified
    // ids use the identity default instance.
    ID compiled_node_id   (ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_element_id(ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_surface_id(ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_line_id   (ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;

    // Part-local geometry construction. All four operations are pre-compile
    // because changing their contents would invalidate dense assembly ids.
    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);
    inline void set_surface(ID id, ID element_id, ID surface_id);
    inline void set_surface(const std::string& elset, ID surface_id);

    // Shared definitions and assignments --------------------------------------
    // Materials, profiles and coordinate systems are independent definition
    // dictionaries and may be extended on either side of compile(). A section
    // added pre-compile belongs to the active Part; a section added post-compile
    // must reference a compiled element set and is assigned immediately.
    void add_coordinate_system(cos::CoordinateSystem::Ptr coordinate_system);
    void add_material(material::Material::Ptr material);
    void add_profile(Profile::Ptr profile);
    void add_section(Section::Ptr section);

    // Loads/supports remain compact registration helpers because both use one
    // uniform data type. Constraints deliberately have no equivalent helper:
    // callers append concrete Connector/Coupling/Tie/Contact/Rbm/Equation
    // objects directly to the matching ModelData collection.
    void add_load     (bc::Load::Ptr load);
    void add_amplitude(bc::Amplitude::Ptr amplitude);
    void add_support  (bc::Support support);

    // adds a point mass to a node with prespecified mass, inertia and spring constants.
    void add_point_mass_feature(const std::string& nset,
                                Precision mass,
                                Vec3 rotary_inertia,
                                Vec3 spring_constants,
                                Vec3 rotary_spring_constants);


    // Element lifecycle and section assignment --------------------------------
    // step_begin()/step_end() own analysis-local caches. They are intentionally
    // separate from compile(), which only creates persistent assembly topology.
    void step_begin();
    void step_end();
    void assign_sections();
    void build_shell_element_normals(Precision equalize_angle_degrees = Precision(20));

    // Material state -----------------------------------------------------------
    [[nodiscard]] Index maximum_material_state_size() const;
    void initialize_material_state(Field& state) const;

    // Global assembly ----------------------------------------------------------
    SystemDofIds build_unconstrained_index_matrix();
    Field build_load_matrix(
        std::vector<std::string> load_sets = {},
        Precision time = 0);
    constraint::ConstraintGroups collect_constraints(
        SystemDofIds& system_dof_ids,
        const std::vector<std::string>& supp_sets = {});
    std::vector<std::pair<bc::Amplitude::Ptr, Field>> build_load_basis(
        std::vector<std::string> load_sets = {});
    SparseMatrix build_stiffness_matrix(
        SystemDofIds& indices,
        const Field* stiffness_scalar = nullptr);
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
    SparseMatrix build_lumped_mass_matrix(
        SystemDofIds& indices);

    // Result recovery ----------------------------------------------------------
    Field compute_stress_state(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl = false);
    Field compute_shell_resultants(Field& displacement);
    Field compute_compliance(Field& displacement);
    Field compute_compliance_angle_derivative(Field& displacement);
    Field compute_volumes();
    Field compute_section_forces(Field& displacement);
    Field compute_shear_flow(Field& displacement);

    // Human-readable model diagnostics. The overview reports semantic topology,
    // compiled assembly data and associated definitions through the logger. The
    // stream operator writes a compact summary only to the supplied stream.
    void print_overview() const;
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);
};

#include "model.ipp"

} // namespace fem::model
