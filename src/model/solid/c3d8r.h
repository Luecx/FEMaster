/**
 * @file c3d8r.h
 * @brief Declares the reduced-integration C3D8 solid with finite-strain physical hourglass stabilization.
 *
 * The continuum response is evaluated at the element center. The twelve
 * non-affine hourglass modes use a parameter-free Belytschko-Bindeman physical
 * stiffness embedded in an objective Total-Lagrangian modal potential.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

#pragma once

#include "c3d8.h"

namespace fem::model {

/**
 * @brief Eight-node reduced-integration hexahedral solid.
 *
 * The continuum contribution uses one material point at the element center.
 * Hourglass stabilization uses the four Flanagan-Belytschko scalar modes after
 * projection against affine displacement fields. Geometry-dependent physical
 * coefficients and the initial isotropic-equivalent shear modulus and Poisson
 * ratio define the Belytschko-Bindeman modal stiffness without a user hourglass
 * coefficient or bulk-modulus singularity.
 *
 * For finite deformation the modal displacement vectors are pulled back by the
 * mean deformation gradient before entering the referential stabilization
 * potential. Residual and tangent are exact first and second derivatives of the
 * same potential, while the undeformed tangent is exactly the physical
 * Belytschko-Bindeman stiffness.
 */
class C3D8R final : public C3D8 {
public:
    static constexpr Index N         = 8;
    static constexpr Dim   D         = 3;
    static constexpr Index n_hg_mode = 4;
    static constexpr Index ndof      = N * D;
    static constexpr Index n_hg_dof  = D * n_hg_mode;

    using GradientMatrix = StaticMatrix<N, D>;
    using HourglassModes = StaticMatrix<N, n_hg_mode>;
    using Matrix24       = StaticMatrix<ndof, ndof>;
    using Vector24       = StaticVector<ndof>;
    using ModalMatrix    = StaticMatrix<n_hg_dof, n_hg_dof>;
    using ModalVector    = StaticVector<n_hg_dof>;
    using ModalJacobian  = StaticMatrix<n_hg_dof, ndof>;

    C3D8R(ID elem_id, const std::array<ID, N>& node_ids);
    ~C3D8R() override = default;

    // Element identification and reduced quadrature
    std::string type_name() const override;
    const math::quadrature::Quadrature& integration_scheme_stiffness() const override;
    RowMatrix stress_strain_nodal_rst() override;

    // Continuum plus hourglass tangent and internal force
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_tangent(Precision* buffer,
                                Field&       ip_stress_state,
                                NodeData&    nodal_forces,
                                const Field& displacement) override;
    void compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) override;

private:
    // Reference physical stabilization data
    HourglassModes  primitive_hourglass_modes();
    GradientMatrix  mean_reference_gradient();
    Mat3            hourglass_reference_frame();
    Mat3            hourglass_geometry_integrals(const Mat3& frame);
    StaticVector<2> hourglass_material_parameters();
    void            build_hourglass_reference_data();
    Matrix24        hourglass_stiffness();

    // Objective finite-strain Total-Lagrangian hourglass response
    void hourglass_response(Vector24& force, Matrix24* tangent);

    GradientMatrix mean_gradient_cache   = GradientMatrix::Zero();
    HourglassModes hourglass_modes_cache = HourglassModes::Zero();
    Mat3           hourglass_frame_cache = Mat3::Identity();
    ModalMatrix    modal_stiffness_cache = ModalMatrix::Zero();
    Matrix24       hourglass_stiffness_cache = Matrix24::Zero();
    bool           hourglass_reference_cached = false;

    // Global force scattering
    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model
