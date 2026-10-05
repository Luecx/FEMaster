/**
 * @file truss.h
 * @brief Declares the three-dimensional two-node truss element.
 *
 * The truss element represents a straight two-node member that carries only
 * axial force. Each node contributes the three translational degrees of freedom
 * of the structural model, while the element kinematics reduce the deformation
 * to a single axial strain measure along the member axis.
 *
 * Linear result recovery uses infinitesimal axial strain and Cauchy stress. The
 * nonlinear formulation follows a Total-Lagrangian description based on the
 * reference length, Green-Lagrange axial strain and second Piola-Kirchhoff
 * stress. The element owns one constitutive material point whose history is
 * addressed through the model-wide committed and trial material-state fields.
 *
 * Mechanical assembly, geometry, field integration and result recovery are
 * implemented in separate translation units around the same compact `T3`
 * interface. Temporary stresses and constitutive tangents remain local to the
 * operation that needs them and are never exposed as solver-level scratch fields.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../core/core.h"
#include "../../material/elasticity.h"
#include "../../material/strain/axial_strain_green_lagrange.h"
#include "../../material/stress/axial_stress_cauchy.h"
#include "../../material/stress/axial_stress_pk2.h"
#include "../../section/section_truss.h"
#include "../element/element_structural.h"

#include <array>
#include <string>

namespace fem {
namespace model {

/**
 * @brief Three-dimensional two-node truss with one axial constitutive point.
 *
 * `T3` is a structural line element with two nodes and three translational
 * degrees of freedom per node. The element has no rotational stiffness and no
 * bending or shear response. Its mechanical behavior is fully determined by the
 * assigned `TrussSection`, the section reference area and the associated
 * material law.
 *
 * The nonlinear kinematics use the reference length `L0` and the length `l0`
 * of the supplied linearization configuration, with `lambda0 = l0 / L0`.
 * Green-Lagrange axial strain follows from this stretch and PK2 stress is
 * evaluated by the constitutive model in the Total-Lagrangian material
 * description. The base stress supplies the geometric part of the complete
 * tangent, while the constitutive tangent supplies its material part.
 *
 * The common mechanical evaluation is expressed about a base state u0. The
 * complete tangent is evaluated at u0, while force at a different state u is
 * obtained by first-order continuation from that base state. The separately
 * requested geometric stiffness is generated only by the linearized stress
 * increment from u0 to u.
 *
 * State-neutral stiffness, perturbation-geometric evaluation and result recovery
 * read committed material history without providing a writable target state.
 * Temporary stresses and constitutive tangents remain local to the operation
 * that requires them.
 *
 * Mechanical evaluation constructs its state explicitly from the reference
 * geometry and the supplied displacement fields. The model's current POSITION
 * field is used only by geometry-based utilities such as volume and distributed
 * field integration.
 */
struct T3 : StructuralElement {
    static constexpr Index N = 2;

    // Element connectivity in global node-id space
    std::array<ID, N> node_ids {};

    T3(ID elem_id, std::array<ID, N> node_ids);
    ~T3() override = default;

    // Recreate only id and connectivity for Instance expansion. Section,
    // offsets, material state and ModelData binding belong to compiled storage.
    ElementPtr copy() const override { return std::make_shared<T3>(elem_id, node_ids); }

    // Element topology and identification. The truss exposes three translational
    // DOFs at each of its two nodes, one material point and no surface topology.
    ElDofs      dofs() const override;
    Dim         dimensions() const override;
    Dim         n_nodes() const override;
    Dim         num_ip() const override;
    const ID*   nodes() const override;
    SurfacePtr  surface(ID surface_id) override;
    std::string type_name() const override;

    // Resolve the assigned truss section, material and elasticity. These access
    // functions validate the compiled element definition before returning the
    // corresponding model objects used by the mechanical operators below.
    TrussSection*         get_section();
    material::MaterialPtr get_material();
    material::Elasticity* get_elasticity();

    // Reference and current geometry. The reference configuration supplies the
    // material coordinate X, reference length L0 and unit axis N0. The current
    // coordinates are used only by geometry-based utilities such as volume and
    // distributed-field integration; mechanical state evaluation constructs
    // its configuration explicitly from X + u0.
    Vec3      node_position_reference(Index local_node) const;
    Vec3      node_position_current(Index local_node) const;
    Precision length_reference() const;
    Precision length_current() const;
    Vec3      direction_reference() const;

    // Current geometric volume A0*l used by generic geometry/integration paths.
    Precision volume() override;

    // Mechanical response through the common state-based interface. A null
    // linearization denotes u0 = 0; a non-null state evaluates the
    // Total-Lagrangian truss tangent at u0. The separate geometric output is
    // generated by the linearized stress increment from u0 to displacement u.
    MapMatrix evaluate(
        Precision*   tangent,
        Precision*   geometric_tangent,
        NodeData*    internal_force,
        const Field* displacement,
        const Field* linearization,
        const Field* thermal_free_strain,
        bool         update_state
    ) override;
    MapMatrix mass(Precision* buffer) override;

    // Natural-coordinate locations used for nodal and integration-point recovery.
    // The axial solution is constant over the element for the supported truss
    // formulation, but both locations are exposed for generic result pipelines.
    RowMatrix stress_strain_nodal_rst() override;
    RowMatrix stress_strain_ip_rst() override;

    // Integrate scalar, vector and tensor fields over the current truss measure
    // A*l. Density scaling is optional and is applied through the assigned
    // material when requested. The nodal vector overload distributes the
    // integrated force equally to the two truss nodes.
    Precision integrate_scalar_field(
        bool               scale_by_density,
        const ScalarField& field
    ) override;
    Vec3 integrate_vector_field(
        bool            scale_by_density,
        const VecField& field
    ) override;
    void integrate_vector_field(
        Field&          node_loads,
        bool            scale_by_density,
        const VecField& field
    ) override;
    Mat3 integrate_tensor_field(
        bool            scale_by_density,
        const TenField& field
    ) override;

    // Thermal expansion from the current nodal temperature state.
    void apply_thermal_expansion_load(
        Field&       node_loads,
        const Field& node_temp
    ) override;
    void apply_thermal_free_strain(
        Field&       thermal_free_strain,
        const Field& node_temp
    ) override;

    // Recover axial strain and physical stress at requested output positions.
    // Linear recovery uses infinitesimal axial strain and Cauchy stress. Finite-
    // strain recovery evaluates Green-Lagrange strain and PK2 stress and pushes
    // the latter forward to axial Cauchy stress for user-facing output.
    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     displacement,
        const RowMatrix& rst,
        int              offset,
        const Field*     linearization,
        const Field*     thermal_free_strain = nullptr
    ) override;
    bool compute_peeq(
        Field& peeq,
        int    offset
    ) override;

    // Element-level scalar compliance contribution based on the current
    // displacement field and the truss stiffness operator.
    void compute_compliance(
        Field& displacement,
        Field& result
    ) override;

    // Recover the constant axial section force at both truss nodes for generic
    // beam-style section-force output. Component zero stores the axial force;
    // remaining output components are cleared by the implementation.
    bool compute_beam_section_forces(
        Field&       section_forces,
        const Field& displacement,
        int          offset,
        const Field* linearization = nullptr
    ) override;
};

} // namespace model
} // namespace fem
