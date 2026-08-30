/**
 * @file surface_shape_matrix.inl
 * @brief Integrates scalar-weighted products of surface shape functions.
 */

#pragma once

#include "surface.h"

namespace fem::model {

template<Index N>
DynamicMatrix Surface<N>::integrate_scalar_shape_matrix(
    const Field&       node_coords,
    const ScalarField& field
) const {
    const auto coordinates = node_coords_global(node_coords);

    // Products N_i N_j need a higher-order rule than ordinary surface loads.
    // Keep this rule local to shape-matrix integration so mechanical pressure,
    // traction and area integration retain their established quadrature.
    const auto& scheme = []() -> const math::quadrature::Quadrature& {
        if constexpr (N == 3) {
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_TRI,
                math::quadrature::ORDER_QUADRATIC
            };
            return q;
        } else if constexpr (N == 4) {
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_QUAD,
                math::quadrature::ORDER_QUADRATIC
            };
            return q;
        } else if constexpr (N == 6) {
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_TRI,
                math::quadrature::ORDER_QUARTIC
            };
            return q;
        } else {
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_QUAD,
                math::quadrature::ORDER_QUARTIC
            };
            return q;
        }
    }();

    StaticMatrix<N, N> matrix = StaticMatrix<N, N>::Zero();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        const Precision r = point.r;
        const Precision s = point.s;

        const StaticMatrix<N, 1> shape = shape_function(r, s);
        const auto jac                   = jacobian(coordinates, r, s);
        const auto position              = interpolate(coordinates, r, s);
        const Precision weighted_area    = jac.col(0).cross(jac.col(1)).norm() * point.w;
        const Precision value            = field(position);

        matrix.noalias() += value * weighted_area * (shape * shape.transpose());
    }

    return matrix;
}

} // namespace fem::model
