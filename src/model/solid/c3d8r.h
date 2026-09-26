/**
 * @file c3d8r.h
 * @brief Declares the reduced-integration C3D8 solid with finite-strain hourglass stabilization.
 *
 * The continuum response is evaluated at the element center. Hourglass control
 * penalizes the non-constant deviatoric Green-Lagrange strain field and the
 * non-uniform volume ratio resolved by a full 2x2x2 geometry sampling. The
 * formulation is objective, uses only the generic material tangent interface
 * and requires no empirical hourglass factor.
 *
 * @author Finn Eggers
 * @date 26.09.2026
 */

#pragma once

#include "c3d8.h"

namespace fem::model {

class C3D8R final : public C3D8 {
public:
    static constexpr Index N    = 8;
    static constexpr Dim   D    = 3;
    static constexpr Index ndof = N * D;

    using Matrix24 = StaticMatrix<ndof, ndof>;
    using Vector24 = StaticVector<ndof>;

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
    Mat6     deviatoric_reference_tangent();
    Matrix24 linear_hourglass_stiffness();
    void     finite_hourglass(
        const Field& displacement,
        Vector24&    local_force,
        Matrix24*    tangent
    );

    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model
