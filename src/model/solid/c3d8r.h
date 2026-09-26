/**
 * @file c3d8r.h
 * @brief Declares the reduced-integration C3D8 solid with physical hourglass stabilization.
 *
 * The continuum response is evaluated at the element center. Four projected
 * Flanagan-Belytschko scalar modes span the hourglass subspace. Stabilization
 * restores the deviatoric stiffness lost by one-point integration inside that
 * subspace instead of assigning all modes one empirical scalar stiffness.
 *
 * @author Finn Eggers
 * @date 26.09.2026
 */

#pragma once

#include "c3d8.h"

namespace fem::model {

/**
 * @brief Eight-node reduced-integration hexahedral solid.
 *
 * The continuum contribution uses one material point at the element center.
 * Hourglass stabilization acts only in the twelve-dimensional subspace formed
 * by four scalar hourglass modes in three translational directions.
 *
 * The stabilization is assembled from the difference between fully integrated
 * and one-point deviatoric reference stiffness, projected onto that subspace.
 * It therefore requires no material-specific modulus and no empirical
 * hourglass coefficient.
 */
class C3D8R final : public C3D8 {
public:
    static constexpr Index N               = 8;
    static constexpr Dim   D               = 3;
    static constexpr Index ndof            = N * D;
    static constexpr Index n_hourglass_dof = 4 * D;

    using GradientMatrix = StaticMatrix<N, D>;
    using HourglassModes = StaticMatrix<N, 4>;
    using HourglassBasis = StaticMatrix<ndof, n_hourglass_dof>;
    using Matrix12       = StaticMatrix<n_hourglass_dof, n_hourglass_dof>;
    using Matrix24       = StaticMatrix<ndof, ndof>;
    using Vector24       = StaticVector<ndof>;

    C3D8R(ID elem_id, const std::array<ID, N>& node_ids);
    ~C3D8R() override = default;

    ElementPtr copy() const override { return std::make_shared<C3D8R>(elem_id, node_ids); }

    std::string type_name() const override;
    const math::quadrature::Quadrature& integration_scheme_stiffness() const override;
    RowMatrix stress_strain_nodal_rst() override;

    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_tangent(
        Precision*   buffer,
        NodeData&    nodal_forces,
        const Field& displacement
    ) override;

protected:
    const RowMatrix& extrapolation_matrix() override {
        static const RowMatrix matrix = math::extrapolate(
            this->stress_strain_ip_rst(), this->node_coords_local(),
            {math::ExtrapolationBasis::F1});
        return matrix;
    }

private:
    HourglassModes primitive_hourglass_modes();
    GradientMatrix mean_reference_gradient();
    HourglassBasis hourglass_basis();
    Mat6           deviatoric_reference_tangent();
    Matrix24       hourglass_stiffness();

    Vector24 local_displacement(const Field& displacement);
    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model
