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
 * Mechanical assembly, geometric quantities, field integration and structural
 * result recovery are implemented directly by `T3`. Temporary stresses and
 * constitutive tangents required during assembly remain local to the element
 * evaluation and are not exposed as solver-level scratch fields.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../core/core.h"
#include "../../material/elasticity.h"
#include "../../material/strain/axial_strain_green_lagrange.h"
#include "../../material/strain/axial_strain_linearized.h"
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
 * The nonlinear kinematics use the reference length `L0`, the current length
 * `l` and the stretch `lambda = l / L0`. Green-Lagrange axial strain is obtained
 * from the stretch, while PK2 stress is evaluated by the constitutive model in
 * the Total-Lagrangian material description. The resulting axial stress drives
 * both the internal force and the geometric part of the tangent matrix. The
 * material tangent from the same constitutive evaluation supplies the material
 * contribution to the consistent nonlinear tangent.
 *
 * The nonlinear equilibrium operator follows one physical constitutive path.
 * `stiffness_tangent()` always evaluates the internal force for the supplied
 * trial configuration. A non-null matrix buffer additionally requests material
 * and geometric tangent assembly; a null buffer performs an internal-force-only
 * evaluation. Both variants read committed material history and write the
 * corresponding persistent trial material state.
 *
 * Linear stiffness, geometric-stiffness evaluation for prestress analyses and
 * result-recovery routines are state-neutral with respect to persistent trial
 * history. They read committed history when required by the constitutive model
 * but pass no target state, so no material-state update is stored. Stress values
 * and constitutive tangents needed only during an element operation remain local
 * to that operation.
 *
 * Reference and current geometry are accessed through the model data bound to
 * the element. The reference axis is used by linearized strain recovery, while
 * the current axis defines the direction of nonlinear axial force and tangent
 * contributions.
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

    // Reference and current geometry. The reference configuration determines
    // L0 and the linearized axial direction; the current configuration determines
    // l, stretch and the force direction used in the nonlinear formulation.
    Vec3      node_position_reference(Index local_node) const;
    Vec3      node_position_current(Index local_node) const;
    Precision length_reference() const;
    Precision length_current() const;
    Precision stretch() const;
    Vec3      direction_reference() const;
    Vec3      direction_current() const;

    // Generic line-element geometry access used by shared structural utilities
    Precision length();
    Vec3      direction();

    // Mechanical operators. `stiffness()` returns the material/linear operator.
    // `stiffness_geom()` derives the stress required for the geometric operator
    // directly from the supplied displacement field without modifying persistent
    // trial history. `stiffness_tangent()` evaluates the nonlinear internal force
    // and assembles the consistent tangent only when `buffer` is non-null.
    Precision volume() override;
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(
        Precision*   buffer,
        const Field& displacement
    ) override;
    MapMatrix stiffness_tangent(
        Precision*   buffer,
        NodeData&    nodal_forces,
        const Field& displacement
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

    // Equivalent nodal loading caused by a prescribed temperature field
    void apply_tload(
        Field&       node_loads,
        const Field& node_temp,
        Precision    ref_temp
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
        bool             use_green_lagrange_nl
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
        int          offset
    ) override;
};

} // namespace model
} // namespace fem
