#pragma once

#include "../../math/interpolate.h"
#include "../element/element_structural.h"
#include "../geometry/surface/surface.h"

#include <array>
#include <functional>

namespace fem {
namespace model {

template<Index N>
struct SolidElement : StructuralElement {
    std::array<ID, N> node_ids {};

protected:
    static constexpr Dim n_strain = 6;
    static constexpr Dim D        = 3;

    friend quadrature::Quadrature;

public:
    SolidElement(ID elem_id, std::array<ID, N> node_ids);
    ~SolidElement() override = default;

    ElDofs dofs                 () override;
    Dim    dimensions           () override;
    Dim    n_nodes              () override;
    Dim    n_integration_points () override;
    ID*    nodes                () override;

    SurfacePtr surface(ID surface_id) override = 0;

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
        const Field& displacement,
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

    SolidSection* get_section();

    bool has_section_orientation () const;
    bool has_topo_orientation    () const;
    bool has_material_orientation() const;

    Vec3 topo_angles              () const;
    Mat3 section_basis            (Precision r, Precision s, Precision t);
    Mat3 topo_basis               () const;
    Mat3 topo_basis_derivative    (Index angle_id) const;
    Mat3 material_basis           (Precision r, Precision s, Precision t);
    Mat3 material_basis_derivative(Precision r, Precision s, Precision t, Index angle_id);

    StaticMatrix<n_strain, n_strain> material_matrix(Precision r, Precision s, Precision t);
    void material_stress_strain(
        Precision                     r,
        Precision                     s,
        Precision                     t,
        const StaticVector<n_strain>& global_strain,
        StaticVector<n_strain>&       out_stress,
        StaticVector<n_strain>&       out_strain
    );

    virtual StaticMatrix<N, 1> shape_function  (Precision r, Precision s, Precision t) = 0;
    virtual StaticMatrix<N, D> shape_derivative(Precision r, Precision s, Precision t) = 0;
    virtual StaticMatrix<N, D> node_coords_local() = 0;
    virtual StaticMatrix<N, D> node_coords_global();

    virtual const quadrature::Quadrature& integration_scheme() = 0;
    virtual const quadrature::Quadrature& integration_scheme_mass();

    template<Dim K>
    StaticVector<K> interpolate(StaticMatrix<N, K> data, Precision r, Precision s, Precision t);

    StaticMatrix<n_strain, D * N> strain_displacement(
        const StaticMatrix<N, D>& shape_der_global
    );
    StaticMatrix<n_strain, D * N> strain_displacements(
        const StaticMatrix<N, D>& node_coords,
        Precision                 r,
        Precision                 s,
        Precision                 t,
        Precision&                det,
        bool                      check_det = true
    );
    StaticMatrix<D, D> jacobian(
        const StaticMatrix<N, D>& node_coords,
        Precision                 r,
        Precision                 s,
        Precision                 t
    );

    template<Dim K>
    StaticMatrix<N, K> nodal_data(const Field& full_data, Index offset = 0, Index stride = 1);

    template<class ElementType>
    static bool test_implementation(bool print = false);
};

} // namespace model
} // namespace fem

#include "element_solid_compute.ipp"
#include "element_solid_load.ipp"
#include "element_solid_.ipp"
#include "element_solid.ipp"
#include "element_solid_test.ipp"
