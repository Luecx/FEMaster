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
 * @date 25.08.2026
 */

#pragma once

#include "../bc/amplitude.h"
#include "../bc/dirichlet/dirichlet.h"
#include "../bc/neumann/neumann.h"
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
 *
 * After compilation the same interface coordinates section assignment,
 * analysis-step element state, material history initialization, global operator
 * assembly and result recovery. Persistent data always remains owned by
 * `ModelData`; `Model` only implements the operations and invariants governing
 * that data.
 *
 * A model is intentionally non-copyable because copied semantic and compiled
 * identifiers could diverge. Move construction transfers the single shared data
 * root without duplicating any finite-element object.
 */
struct Model {
    // Root-level unqualified topology is stored as an implicit Part and embedded
    // through an identity Instance. This lets all semantic entity references use
    // one consistent namespace rule without a separate Assembly object.
    static constexpr const char* DEFAULT_PART_NAME     = "__DEFAULT_PART__";
    static constexpr const char* DEFAULT_INSTANCE_NAME = "__DEFAULT_INSTANCE__";

    // Single source of truth for semantic definitions, reusable topology and
    // compiled solver-facing storage. All Model operations act directly on this
    // shared root; no mirrored compile state or secondary entity container is
    // maintained by the behavioral interface.
    ModelDataPtr _data;

    // Construction and ownership. The default constructor creates the implicit
    // root Part and identity Instance used for unqualified input entities.
    // Copying is forbidden, while move construction transfers the data root and
    // every object reachable from it.
    Model();
    Model(const Model&)            = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) noexcept        = default;

    // Semantic topology ownership. Parts contain reusable local geometry and
    // regions; instances embed those definitions through rigid placements.
    // Both operations are valid only before the one-way compilation boundary.
    void add_part(const std::string& name);
    void add_instance(const std::string& name,
                      const std::string& part,
                      Vec3 translation = Vec3::Zero(),
                      Mat3 rotation = Mat3::Identity());

    // One-way semantic-to-assembly transition. Compilation creates deterministic
    // dense identifiers, transforms nodal geometry and section orientations,
    // rewires copied element/surface/line connectivity, materializes inherited
    // regions and initializes element-nodal, integration-point and material-point
    // enumeration. The operation may be called exactly once.
    void compile();

    // Resolution of semantic local identifiers after compilation. The typed
    // overloads address a specified Instance map. String overloads accept
    // either `ID` or `INSTANCE.ID` and use the implicit identity Instance for
    // unqualified references.
    ID compiled_node_id   (ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_element_id(ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_surface_id(ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;
    ID compiled_line_id   (ID id, const std::string& instance = DEFAULT_INSTANCE_NAME) const;

    ID compiled_node_id   (const std::string& reference) const;
    ID compiled_element_id(const std::string& reference) const;
    ID compiled_surface_id(const std::string& reference) const;

    // Named-region population in explicit semantic identifier spaces. Part
    // operations use the active Part and retain sparse local identifiers;
    // assembly operations use compiled regions and map scalar `ID` or
    // `INSTANCE.ID` references into dense identifiers.
    void add_nodes_to_part_set         (const std::string& set, const std::string& source);
    void add_nodes_to_assembly_set     (const std::string& set, const std::string& source);
    void add_elements_to_part_set      (const std::string& set, const std::string& source);
    void add_elements_to_assembly_set  (const std::string& set, const std::string& source);
    void add_surfaces_to_part_set      (const std::string& set, const std::string& source);
    void add_surfaces_to_assembly_set  (const std::string& set, const std::string& source);

    // Part-local topology construction in the currently active semantic Part.
    // Nodes and element connectivity retain their input identifiers, while
    // boundary extraction creates a local surface or line from an element side.
    // These operations are pre-compile because later changes would invalidate
    // dense assembly identifiers and inherited regions.
    inline void set_node(ID id, Precision x, Precision y, Precision z = 0);
    template<typename T, typename... Args>
    inline void set_element(ID id, Args&&... args);
    inline void set_surface(ID id, ID element_id, ID surface_id);
    inline void set_surface(const std::string& elset, ID surface_id);

    // Shared definitions and section assignments. Materials, profiles and
    // coordinate systems are independent of identifier flattening and may be
    // registered on either side of compile(). A pre-compile section belongs to
    // the active Part and is copied per Instance; a post-compile section must
    // already reference a compiled element region and is assigned immediately.
    void add_coordinate_system(cos::CoordinateSystem::Ptr coordinate_system);
    void add_material(material::Material::Ptr material);
    void add_profile(Profile::Ptr profile);
    void add_section(Section::Ptr section);

    // Boundary-condition resources and collector entries. Loads and supports
    // address compiled regions and therefore require an active collector after
    // compilation. Amplitudes are shared named definitions and remain independent
    // of the topology transition. Constraints deliberately have no registration
    // helper; callers append their concrete type directly to ModelData.
    void add_load     (bc::Neumann::Ptr load);
    void add_amplitude(bc::Amplitude::Ptr amplitude);
    void add_support  (bc::Dirichlet::Ptr support);

    // Compiled element preparation and analysis lifecycle. Section assignment
    // binds compiled elements to their section definitions, and shell-normal
    // construction prepares element-nodal reference geometry. step_begin() and
    // step_end() manage analysis-local element caches independently of the
    // persistent topology created by compile().
    void step_begin();
    void step_end();
    void assign_sections();
    void build_shell_element_normals(Precision equalize_angle_degrees = Precision(20));

    // Material-point history allocation and initialization. The maximum width is
    // derived from the elastic laws assigned to compiled elements; initialization
    // then prepares every enumerated material point in a compatible field.
    [[nodiscard]] Index maximum_material_state_size() const;
    void initialize_material_state(Field& state) const;

    // Global system construction. These operations enumerate active DOFs,
    // collect prescribed loads and constraints, and assemble structural tangent,
    // geometric, mass and internal-force contributions in global coordinates.
    // Optional stiffness-scaling fields act per element.
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
    SparseMatrix build_lumped_mass_matrix(SystemDofIds& indices);

    // Result recovery from compiled structural elements. Returned fields use the
    // domain appropriate to each quantity, including integration-point stress,
    // element-nodal stress/strain, shell resultants and element-level compliance,
    // volume, section-force and shear-flow measures.
    Field compute_stress_state(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl = false);
    std::tuple<Field, Field> compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl = false);
    Field compute_shell_resultants(Field& displacement);
    Field compute_compliance(Field& displacement);
    Field compute_compliance_angle_derivative(Field& displacement);
    Field compute_volumes();
    Field compute_section_forces(Field& displacement);
    Field compute_shear_flow(Field& displacement);

    // Human-readable diagnostics. The overview reports semantic topology,
    // compiled assembly data and associated definitions through the hierarchical
    // logger. The stream operator writes a compact, side-effect-free summary only
    // to the supplied stream.
    void print_overview() const;
    friend std::ostream& operator<<(std::ostream& ostream, const Model& model);
};

#include "model.ipp"

} // namespace fem::model
