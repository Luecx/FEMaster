/**
 * @file surface_shape_matrix.inl
 * @brief Implements scalar-weighted surface shape-function product integration.
 *
 * Mixed scalar boundary conditions require local operators of the form
 *
 *     M_s = integral_Gamma a(x) N N^T dGamma.
 *
 * Products of interpolation functions have higher polynomial order than the
 * ordinary surface-load integrands. This implementation therefore selects a
 * dedicated quadrature rule based on interpolation order instead of changing the
 * established quadrature used by pressure, traction and area integration.
 *
 * @see Surface
 * @see SurfaceInterface::integrate_scalar_shape_matrix
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include "surface.h"

namespace fem::model {

/**
 * Integrates a scalar-weighted consistent surface shape matrix.
 *
 * Linear three- and four-node surfaces use quadratic quadrature because
 * `N_i N_j` is quadratic. Quadratic six- and eight-node surfaces use quartic
 * quadrature. At every point the physical area measure is
 *
 *     dGamma = ||x_,r cross x_,s|| w_q,
 *
 * and the local contribution is
 *
 *     a(x_q) N(x_q) N(x_q)^T dGamma.
 *
 * @param node_coords Global nodal coordinates defining the physical surface.
 * @param field Scalar coefficient evaluated at global quadrature positions.
 * @return Dense local matrix in surface-connectivity ordering.
 */
template<Index N>
DynamicMatrix Surface<N>::integrate_scalar_shape_matrix(
    const Field&       node_coords,
    const ScalarField& field
) const {
    // Gather physical geometry once for all quadrature points
    const auto coordinates = node_coords_global(node_coords);

    // Use a dedicated rule for the higher polynomial order of N_i N_j
    const auto& scheme = []() -> const math::quadrature::Quadrature& {
        if constexpr (N == 3) {
            static const math::quadrature::Quadrature quadrature{
                math::quadrature::DOMAIN_ISO_TRI,
                math::quadrature::ORDER_QUADRATIC
            };
            return quadrature;
        } else if constexpr (N == 4) {
            static const math::quadrature::Quadrature quadrature{
                math::quadrature::DOMAIN_ISO_QUAD,
                math::quadrature::ORDER_QUADRATIC
            };
            return quadrature;
        } else if constexpr (N == 6) {
            static const math::quadrature::Quadrature quadrature{
                math::quadrature::DOMAIN_ISO_TRI,
                math::quadrature::ORDER_QUARTIC
            };
            return quadrature;
        } else {
            static const math::quadrature::Quadrature quadrature{
                math::quadrature::DOMAIN_ISO_QUAD,
                math::quadrature::ORDER_QUARTIC
            };
            return quadrature;
        }
    }();

    StaticMatrix<N, N> matrix = StaticMatrix<N, N>::Zero();

    // Integrate a(x) N N^T over the complete physical surface
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        const StaticMatrix<N, 1> shape    = shape_function(point.r, point.s);
        const auto               jac      = jacobian(coordinates, point.r, point.s);
        const auto               position = interpolate(coordinates, point.r, point.s);

        const Precision weighted_area =
            jac.col(0).cross(jac.col(1)).norm() * point.w;
        const Precision value = field(position);

        matrix.noalias() += value * weighted_area * (shape * shape.transpose());
    }

    return matrix;
}

} // namespace fem::model
