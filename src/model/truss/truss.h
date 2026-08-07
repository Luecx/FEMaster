/**
 * @file truss.h
 * @brief Declares the three-dimensional two-node truss element.
 *
 * The truss supports linearized axial output and a Total-Lagrangian nonlinear
 * formulation based on Green-Lagrange strain and PK2 stress. Its single
 * constitutive material point is addressed directly through the model material
 * state field.
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
 * The nonlinear tangent evaluates axial PK2 stress and its material tangent in
 * one constitutive call. The same stress drives geometric stiffness and internal
 * force, preventing an in-place history state from advancing twice within one
 * solver evaluation.
 */
struct T3 : StructuralElement {
    static constexpr Index N = 2;

    std::array<ID, N> node_ids {};

    T3(ID elem_id, std::array<ID, N> node_ids);
    ~T3() override = default;

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

    // Element matrices and nonlinear tangent. The complete tangent evaluates the
    // single material-point state exactly once and reuses its PK2 stress for the
    // geometric block, stored IP stress and matching internal force.
    Precision volume() override;
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;
    MapMatrix stiffness_tangent(Precision* buffer,
                                Field&       ip_stress_state,
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
    // nonlinear path stores PK2 stress for Total-Lagrangian assembly but pushes
    // it forward to axial Cauchy stress for user output. All constitutive calls
    // use the element's single globally enumerated material-state row.
    void compute_stress_strain(Field*           strain,
                               Field*           stress,
                               const Field&     displacement,
                               const RowMatrix& rst,
                               int              offset,
                               bool             use_green_lagrange_nl) override;
    void compute_stress_state(Field&       stress_state,
                              const Field& displacement,
                              int          offset,
                              bool         use_green_lagrange_nl) override;
    void compute_internal_force_nonlinear(Field&       node_forces,
                                          const Field& ip_stress) override;
    void compute_compliance(Field& displacement, Field& result) override;
    bool compute_beam_section_forces(Field&       section_forces,
                                     const Field& displacement,
                                     int          offset) override;
};

} // namespace model
} // namespace fem
