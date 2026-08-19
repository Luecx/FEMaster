/**
 * @file c3d8r.h
 * @brief Declares the reduced-integration C3D8 solid with hourglass stabilization.
 *
 * The continuum response is evaluated at the element center. Four projected
 * Flanagan-Belytschko hourglass modes add a reference-configuration
 * stabilization tangent and matching internal force.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "c3d8.h"

namespace fem::model {

/**
 * @brief Eight-node reduced-integration hexahedral solid.
 *
 * The continuum contribution uses one material point at the element center.
 * Hourglass stabilization acts independently in the three translational
 * directions through
 *
 *     K_hg = k_hg kron(G G^T, I_3),
 *
 * with matching force `f_hg = K_hg u_e`. Auxiliary zero-strain evaluation of
 * the hourglass modulus is state-neutral; the physical center-point material
 * state is advanced only by the continuum constitutive update.
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

    static constexpr Precision default_hourglass_coefficient = Precision(3.0);

    C3D8R(ID elem_id, const std::array<ID, N>& node_ids);
    ~C3D8R() override = default;

    // Preserve the reduced-integration concrete type while copying only the
    // persistent topology. Hourglass/runtime state is rebuilt by the normal
    // compiled-element lifecycle.
    ElementPtr copy() const override { return std::make_shared<C3D8R>(elem_id, node_ids); }

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
    // Hourglass modes, reference gradients and material scaling
    HourglassModes primitive_hourglass_modes();
    GradientMatrix mean_reference_gradient(Precision& reference_volume);
    Precision      hourglass_material_scale();
    Matrix24       hourglass_stiffness();

    // Element-local displacement and global force scattering
    Vector24 local_displacement();
    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model