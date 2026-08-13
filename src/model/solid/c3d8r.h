/**
 * @file c3d8r.h
 * @brief Declares the reduced-integration C3D8 solid with physical hourglass stabilization.
 *
 * The continuum response is evaluated at the element center. The twelve
 * non-affine hourglass modes are stabilized by a parameter-free assumed-strain
 * stiffness following the Belytschko-Bindeman physical stabilization form.
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
 * ratio define the Belytschko-Bindeman assumed-strain stabilization matrix
 * without a user hourglass coefficient or bulk-modulus singularity.
 *
 * The stabilization is formed in an element-fixed referential frame and rotated
 * back to global translational DOFs. The resulting matrix is constant in the
 * reference configuration and supplies the matching force `f_hg = K_hg u_e`.
 */
class C3D8R final : public C3D8 {
public:
    static constexpr Index N    = 8;
    static constexpr Dim   D    = 3;
    static constexpr Index ndof = N * D;

    using GradientMatrix = StaticMatrix<N, D>;
    using HourglassModes = StaticMatrix<N, 4>;
    using Matrix24       = StaticMatrix<ndof, ndof>;
    using Vector24       = StaticVector<ndof>;

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
    // Hourglass basis, referential frame, geometry factors and material parameters
    HourglassModes  primitive_hourglass_modes();
    GradientMatrix  mean_reference_gradient();
    Mat3            hourglass_reference_frame();
    Mat3            hourglass_geometry_integrals(const Mat3& frame);
    StaticVector<2> hourglass_material_parameters();
    Matrix24        hourglass_stiffness();

    Matrix24 hourglass_stiffness_cache = Matrix24::Zero();
    bool     hourglass_stiffness_cached = false;

    // Element-local displacement and global force scattering
    Vector24 local_displacement();
    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model
