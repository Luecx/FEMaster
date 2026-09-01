/**
 * @file truss.h
 * @brief Declares the three-dimensional two-node truss element.
 *
 * The truss supports linearized axial output and a Total-Lagrangian nonlinear
 * formulation based on Green-Lagrange strain and PK2 stress. Its single
 * constitutive material point is addressed directly through the model material
 * input/output state fields.
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
 * @brief Three-dimensional two-node truss with one axial material point.
 *
 * The nonlinear formulation evaluates Green-Lagrange axial strain and PK2
 * stress in the Total-Lagrangian setting. The same constitutive response drives
 * the internal force, material tangent and geometric tangent so all nonlinear
 * contributions correspond to one physical trial state.
 *
 * Material history is stored in the model-wide committed and trial material
 * state fields. Nonlinear equilibrium evaluation reads the committed state and
 * writes the trial state. Linear stiffness and geometric-stiffness evaluations
 * use local constitutive scratch storage and therefore do not advance persistent
 * material history.
 */
struct T3 : StructuralElement {
    static constexpr Index N = 2;

    std::array<ID, N> node_ids {};

    T3(ID elem_id, std::array<ID, N> node_ids);
    ~T3() override = default;

    // Recreate only id and connectivity for Instance expansion. Section,
    // offsets, material state and ModelData binding belong to compiled storage.
    ElementPtr copy() const override { return std::make_shared<T3>(elem_id, node_ids); }

    // Element topology and identification
    ElDofs    dofs() const override;
    Dim       dimensions() const override;
    Dim       n_nodes() const override;
    Dim       num_ip() const override;
    const ID* nodes() const override;
    SurfacePtr surface(ID surface_id) override;
    std::string type_name() const override;

    // Resolve the assigned truss section, material and elasticity, then expose
    // reference/current length, stretch and unit directions. Reference geometry
    // defines Green-Lagrange strain; current direction defines force orientation.
    TrussSection*         get_section();
    material::MaterialPtr get_material();
    material::Elasticity* get_elasticity();

    Vec3      node_position_reference(Index local_node) const;
    Vec3      node_position_current(Index local_node) const;
    Precision length_reference() const;
    Precision length_current() const;
    Precision stretch() const;
    Vec3      direction_reference() const;
    Vec3      direction_current() const;

    Precision length();
    Vec3      direction();

    // Mechanical element operators. The geometric stiffness obtains its PK2
    // stress locally from the supplied displacement state. The nonlinear tangent
    // always assembles the matching internal force and skips tangent assembly
    // when the caller supplies a null matrix buffer.
    Precision volume() override;
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& displacement) override;
    MapMatrix stiffness_tangent(Precision*   buffer,
                                NodeData&    nodal_forces,
                                const Field& displacement) override;
    MapMatrix mass(Precision* buffer) override;

    // Stress/strain recovery coordinates
    RowMatrix stress_strain_nodal_rst() override;
    RowMatrix stress_strain_ip_rst() override;

    // Distributed field integration
    Precision integrate_scalar_field(bool               scale_by_density,
                                     const ScalarField& field) override;
    Vec3      integrate_vector_field(bool            scale_by_density,
                                     const VecField& field) override;
    void      integrate_vector_field(Field&          node_loads,
                                     bool            scale_by_density,
                                     const VecField& field) override;
    Mat3      integrate_tensor_field(bool            scale_by_density,
                                     const TenField& field) override;

    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override;

    // Recover strain and physical stress at requested output positions. The
    // nonlinear output path evaluates PK2 stress and pushes it forward to axial
    // Cauchy stress for user output. Solver-internal stress quantities remain
    // local to the corresponding mechanical operator.
    void compute_stress_strain(Field*           strain,
                               Field*           stress,
                               const Field&     displacement,
                               const RowMatrix& rst,
                               int              offset,
                               bool             use_green_lagrange_nl) override;
    void compute_compliance(Field& displacement, Field& result) override;
    bool compute_beam_section_forces(Field&       section_forces,
                                     const Field& displacement,
                                     int          offset) override;
};

} // namespace model
} // namespace fem
