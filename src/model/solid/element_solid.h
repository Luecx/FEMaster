/**
 * @file element_solid.h
 * @brief Declares the common three-dimensional solid element formulation.
 *
 * `SolidElement` provides reference/current geometry, material evaluation,
 * strain-displacement operators and the common linear, prestressed and
 * Total-Lagrangian nonlinear assembly used by the concrete C3D solid topologies.
 *
 * Every constitutive integration point is associated with one globally
 * enumerated material-point row. State-neutral operators read the committed row
 * without producing a persistent target state. The physical nonlinear tangent
 * evaluation reads committed history and writes the corresponding trial row.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../math/interpolate.h"
#include "../../material/strain/volume_strain_green_lagrange.h"
#include "../../material/strain/volume_strain_linearized.h"
#include "../../material/stress/volume_stress_cauchy.h"
#include "../../material/stress/volume_stress_pk2.h"
#include "../element/element_structural.h"
#include "../element/element_thermal.h"
#include "../geometry/surface/surface.h"

#include <array>
#include <functional>

namespace fem::model {

/**
 * @brief Common base for three-dimensional continuum elements.
 *
 * Concrete solid elements provide topology-specific interpolation and
 * quadrature. The base implements geometry transformations, constitutive
 * evaluation and the common structural operators. `stiffness()` is the linear
 * small-strain operator in the reference configuration. `stiffness_geom()`
 * derives prestress directly from the supplied displacement field without
 * modifying persistent material history. `stiffness_tangent()` performs the
 * physical Total-Lagrangian trial evaluation, assembling internal force in all
 * cases and the complete tangent only when matrix storage is requested.
 *
 * Stress is kept local to the element evaluation. No global integration-point
 * stress scratch field is required for geometric stiffness or nonlinear force
 * assembly. Persistent constitutive history consists only of committed and trial
 * material-state rows.
 */
template<Index N>
struct SolidElement : StructuralElement, ThermalElement{
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
    bool      is_solid() const override { return true; }

    SurfacePtr surface(ID surface_id) override = 0;

    // Element geometry and interpolation. Concrete topologies define natural
    // node coordinates, shape functions and quadrature; the base gathers
    // reference/current global coordinates from ModelData.
    virtual StaticMatrix<N, D> node_coords_local() = 0;
    virtual StaticMatrix<N, D> node_coords_reference();
    virtual StaticMatrix<N, D> node_coords_current();

    virtual StaticMatrix<N, 1> shape_function(Precision r, Precision s, Precision t) = 0;
    virtual StaticMatrix<N, D> shape_derivative(Precision r, Precision s, Precision t) = 0;

    virtual const math::quadrature::Quadrature& integration_scheme() const = 0;
    virtual const math::quadrature::Quadrature& integration_scheme_stiffness() const {
        return integration_scheme();
    }

    // Resolve the assigned solid section and optional element-level material
    // rotation/stiffness scale used by every constitutive evaluation.
    SolidSection* get_section();

    Mat3      additional_material_rotation() const;
    Precision element_stiffness_scale() const;
    Vec3      material_position_reference(Precision r, Precision s, Precision t);

    // Evaluate the zero-Green-Lagrange material tangent at one natural point.
    // A null new_state requests a state-neutral auxiliary constitutive query.
    Mat6 material_tangent_reference(Precision r, Precision s, Precision t,
                                    const Precision* old_state, Precision* new_state);

    // Evaluate linearized Cauchy stress in global coordinates. The material
    // tangent is optional and topology stiffness scaling is applied to every
    // requested output.
    void evaluate_material(Precision                     r,
                           Precision                     s,
                           Precision                     t,
                           const VolumeStrainLinearized& global_strain,
                           const Precision*              old_state,
                           Precision*                    new_state,
                           VolumeStressCauchy&           global_stress,
                           Mat6*                         global_tangent = nullptr);

    // Evaluate Total-Lagrangian PK2 stress in global reference coordinates. A
    // null tangent propagates through section and material evaluation so
    // residual-only assembly does not construct dS/dE.
    void evaluate_material(Precision                        r,
                           Precision                        s,
                           Precision                        t,
                           const VolumeStrainGreenLagrange& global_strain,
                           const Precision*                 old_state,
                           Precision*                       new_state,
                           VolumeStressPK2&                 global_stress,
                           Mat6*                            global_tangent = nullptr);

    // Interpolation and geometry transformations
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

    // Mechanical operators. Linear and prestress operators are state-neutral.
    // The nonlinear tangent performs one physical constitutive trial update per
    // stiffness quadrature point and reuses the local PK2 stress for residual and
    // geometric stiffness assembly.
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

    // Thermal operators
    MapMatrix conductivity(Precision* buffer) override;
    MapMatrix capacity(Precision* buffer) override;

    Precision volume() override;

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

    // Stress/strain recovery is state-neutral. Requested output coordinates may
    // reuse committed history from the nearest constitutive integration point;
    // nonlinear PK2 stress is pushed forward to Cauchy stress for user output.
    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     displacement,
        const RowMatrix& rst,
        int              offset,
        bool             use_green_lagrange_nl) override;
    void compute_compliance(
        Field& displacement,
        Field& result) override;
    void compute_compliance_angle_derivative(
        Field& displacement,
        Field& result) override;
    void compute_heat_flux(
        Field& heat_flux,
        const Field& temperature) override;

    template<class ElementType>
    static bool test_implementation(bool print = false);
};

} // namespace fem::model

#include "element_solid_compute.ipp"
#include "element_solid_load.ipp"
#include "element_solid_.ipp"
#include "element_solid.ipp"
#include "element_solid_test.ipp"