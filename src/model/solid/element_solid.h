#pragma once

#include "../../math/interpolate.h"
#include "../../material/strain/volume_strain_green_lagrange.h"
#include "../../material/strain/volume_strain_linearized.h"
#include "../../material/stress/volume_stress_cauchy.h"
#include "../../material/stress/volume_stress_pk2.h"
#include "../element/element_structural.h"
#include "../geometry/surface/surface.h"

#include <array>
#include <functional>

namespace fem::model {

template<Index N>
struct SolidElement : StructuralElement {
    static constexpr Dim D        = 3;
    static constexpr Dim n_strain = 6;

    std::array<ID, N> node_ids {};

protected:
    friend math::quadrature::Quadrature;

public:
    SolidElement(ID elem_id, std::array<ID, N> node_ids)
        : StructuralElement(elem_id)
        , node_ids(node_ids) {}

    ~SolidElement() override = default;

    ElDofs    dofs() const override;
    Dim       dimensions() const override;
    Dim       n_nodes() const override;
    Dim       num_ip() const override;
    const ID* nodes() const override;
    bool      is_solid () const override {return true;}

    SurfacePtr surface(ID surface_id) override = 0;

    virtual StaticMatrix<N, D> node_coords_local    () = 0;
    virtual StaticMatrix<N, D> node_coords_reference();
    virtual StaticMatrix<N, D> node_coords_current  ();

    virtual StaticMatrix<N, 1> shape_function  (Precision r, Precision s, Precision t) = 0;
    virtual StaticMatrix<N, D> shape_derivative(Precision r, Precision s, Precision t) = 0;

    virtual const math::quadrature::Quadrature& integration_scheme() const = 0;
    virtual const math::quadrature::Quadrature& integration_scheme_mass() const {
        return integration_scheme();
    }

    SolidSection* get_section();

    Mat3      additional_material_rotation() const;
    Precision element_stiffness_scale      () const;
    Vec3      material_position_reference  (Precision r, Precision s, Precision t);

    Mat6 material_tangent_reference(Precision r, Precision s, Precision t);

    void evaluate_material(Precision                     r,
                           Precision                     s,
                           Precision                     t,
                           const VolumeStrainLinearized& global_strain,
                           VolumeStressCauchy&           global_stress,
                           Mat6&                         global_tangent);

    void evaluate_material(Precision                        r,
                           Precision                        s,
                           Precision                        t,
                           const VolumeStrainGreenLagrange& global_strain,
                           VolumeStressPK2&                 global_stress,
                           Mat6&                            global_tangent);

    template<Dim K>
    StaticVector<K> interpolate(StaticMatrix<N, K> data,
                                Precision          r,
                                Precision          s,
                                Precision          t);

    template<Dim K>
    StaticMatrix<N, K> nodal_data(const Field& full_data,
                                  Index        offset = 0,
                                  Index        stride = 1);

    StaticMatrix<D, D> jacobian(const StaticMatrix<N, D>& node_coords,
                                Precision                 r,
                                Precision                 s,
                                Precision                 t);

    Mat3 deformation_gradient(const StaticMatrix<N, D>& reference_coords,
                              const StaticMatrix<N, D>& current_coords,
                              Precision                 r,
                              Precision                 s,
                              Precision                 t);

    StaticMatrix<n_strain, D * N> strain_displacement(
        const StaticMatrix<N, D>& shape_der_global
    );
    StaticMatrix<N, D> shape_derivatives_reference(
        const StaticMatrix<N, D>& reference_coords,
        Precision                 r,
        Precision                 s,
        Precision                 t,
        Precision&                det,
        bool                      check_det = true
    );
    StaticMatrix<n_strain, D * N> green_lagrange_strain_displacement(
        const StaticMatrix<N, D>& dN_dX,
        const Mat3&               F
    );
    StaticMatrix<n_strain, D * N> strain_displacements(
        const StaticMatrix<N, D>& node_coords,
        Precision                 r,
        Precision                 s,
        Precision                 t,
        Precision&                det,
        bool                      check_det = true
    );

    MapMatrix stiffness     (Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;
    MapMatrix mass          (Precision* buffer) override;

    Precision volume() override;

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
        const Field& ip_stress
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

    template<class ElementType>
    static bool test_implementation(bool print = false);
};

} // namespace fem::model

#include "element_solid_compute.ipp"
#include "element_solid_load.ipp"
#include "element_solid_.ipp"
#include "element_solid.ipp"
#include "element_solid_test.ipp"
