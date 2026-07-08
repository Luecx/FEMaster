#pragma once

#include "../../core/core.h"
#include "../../material/isotropic_elasticity.h"
#include "../../section/section_truss.h"
#include "../element/element_structural.h"

#include <array>
#include <string>

namespace fem {
namespace model {

struct T3 : StructuralElement {
    static constexpr Index N = 2;

    std::array<ID, N> node_ids {};

    T3(ID elem_id, std::array<ID, N> node_ids);
    ~T3() override = default;

    ElDofs    dofs() const override;
    Dim       dimensions() const override;
    Dim       n_nodes() const override;
    Dim       n_integration_points() const override;
    const ID* nodes() const override;

    SurfacePtr surface(ID surface_id) override;

    std::string type_name() const override;

    TrussSection*                  get_section();
    material::MaterialPtr          get_material();
    material::IsotropicElasticity* get_elasticity();

    Vec3      node_position_reference(Index local_node) const;
    Vec3      node_position_current  (Index local_node) const;
    Precision length_reference      () const;
    Precision length_current        () const;
    Precision stretch               () const;
    Vec3      direction_reference   () const;
    Vec3      direction_current     () const;

    Precision material_tangent_reference();
    Precision material_tangent_spatial  ();

    Precision length();
    Vec3      direction();

    Precision volume        () override;
    MapMatrix stiffness     (Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;
    MapMatrix mass          (Precision* buffer) override;

    RowMatrix stress_strain_nodal_rst() override;
    RowMatrix stress_strain_ip_rst   () override;

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

    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     displacement,
        const RowMatrix& rst,
        int              offset,
        bool             use_green_lagrange_nl
    ) override;
    void compute_stress_state(
        Field&       stress_state,
        const Field& displacement,
        int          offset,
        bool         use_green_lagrange_nl
    ) override;
    void compute_internal_force_nonlinear(
        Field&       node_forces,
        const Field& ip_stress,
        int          ip_offset
    ) override;
    void compute_compliance(
        Field& displacement,
        Field& result
    ) override;
    void compute_compliance_angle_derivative(
        Field& displacement,
        Field& result
    ) override;
    bool compute_shear_flow(
        Field&       shear_flow,
        const Field& displacement,
        int          offset
    ) override;
    bool compute_beam_section_forces(
        Field&       section_forces,
        const Field& displacement,
        int          offset
    ) override;
    bool compute_shell_section_forces(
        Field&       section_forces,
        Field&       contribution_count,
        const Field& displacement
    ) override;
};

} // namespace model
} // namespace fem
