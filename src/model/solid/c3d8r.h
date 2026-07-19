/**
 * @file c3d8r.h
 * @brief Eight-node reduced-integration hexahedral solid with
 *        Flanagan-Belytschko hourglass stabilization.
 */

#pragma once

#include "c3d8.h"

namespace fem::model {

/**
 * @brief Eight-node reduced-integration hexahedral solid.
 *
 * The continuum contribution is evaluated at the element center. Four scalar
 * hourglass modes are projected against the affine displacement field and
 * stabilized independently in all three translational directions.
 *
 * The hourglass formulation is fixed in the reference configuration:
 *
 *     K_hg = k_hg kron(G G^T, I_3)
 *
 * and
 *
 *     f_hg = K_hg u_e.
 *
 * The scalar stiffness uses the initial constitutive shear scale
 *
 *     G_eff = (C44 + C55 + C66) / 3.
 *
 * Consequently, the hourglass residual and tangent are exactly consistent
 * and remain constant throughout a nonlinear analysis.
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

    /**
     * @brief Dimensionless coefficient applied to the hourglass stiffness.
     */
    static constexpr Precision default_hourglass_coefficient = Precision(3.0);

    C3D8R(ID elem_id, const std::array<ID, N>& node_ids);

    ~C3D8R() override = default;

    std::string type_name() const override;

    /**
     * @brief One-point integration at the element center with weight eight.
     */
    const math::quadrature::Quadrature& integration_scheme() const override;

    /**
     * @brief Evaluates all extrapolated nodal results at the element center.
     *
     * Nonlinear stress output is converted from PK2 to Cauchy stress by the
     * inherited SolidElement result path.
     */
    RowMatrix stress_strain_nodal_rst() override;

    /**
     * @brief Returns the one-point continuum tangent plus hourglass tangent.
     */
    MapMatrix stiffness(Precision* buffer) override;

    /**
     * @brief Assembles one-point continuum and hourglass internal forces.
     */
    void compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) override;

private:
    /**
     * @brief Returns the four primitive scalar hourglass modes.
     */
    HourglassModes primitive_hourglass_modes();

    /**
     * @brief Computes the volume-averaged reference shape gradients.
     *
     * The averaging uses the full 2x2x2 C3D8 integration scheme and is
     * independent of the one-point continuum integration rule.
     */
    GradientMatrix mean_reference_gradient(Precision& reference_volume);

    /**
     * @brief Returns the initial mean constitutive shear stiffness.
     *
     * The material is evaluated at zero Green-Lagrange strain, corresponding
     * to the undeformed state F = I.
     */
    Precision hourglass_material_scale();

    /**
     * @brief Builds the constant 24x24 reference hourglass tangent.
     */
    Matrix24 hourglass_stiffness();

    /**
     * @brief Extracts the local displacement vector in node-major order.
     */
    Vector24 local_displacement();

    /**
     * @brief Adds a local 24-vector to a global nodal field.
     */
    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model
