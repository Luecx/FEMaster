/**
 * @file c3d8r.h
 * @brief Declares the reduced-integration C3D8 solid with projected hourglass stabilization.
 *
 * The continuum response is evaluated at the element center. The twelve
 * non-affine hourglass modes are stabilized by the fully integrated deviatoric
 * reference stiffness projected onto their subspace.
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
 * Hourglass stabilization is parameter-free and acts only on the non-affine
 * displacement subspace through
 *
 *     K_hg = P_hg K_dev^(8GP) P_hg,
 *
 * where `P_hg` is the orthogonal projector onto the twelve projected hourglass
 * modes and `K_dev^(8GP)` is the fully integrated deviatoric reference stiffness.
 * The resulting matrix is constant in the reference configuration and supplies
 * the matching force `f_hg = K_hg u_e`.
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
    // Hourglass basis and parameter-free projected reference stiffness
    HourglassModes primitive_hourglass_modes();
    GradientMatrix mean_reference_gradient();
    Matrix24       hourglass_stiffness();

    Matrix24 hourglass_stiffness_cache = Matrix24::Zero();
    bool     hourglass_stiffness_cached = false;

    // Element-local displacement and global force scattering
    Vector24 local_displacement();
    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model
