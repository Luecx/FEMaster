/**
 * @file surface_shape_matrix.inl
 * @brief Implements scalar-weighted products of surface shape functions.
 *
 * The routines in this file integrate matrices of the form
 *
 * \f[
 *     M_s = \int_\Gamma a(\mathbf{x})\,N N^T\,\mathrm{d}\Gamma,
 * \f]
 *
 * on the physical finite-element surface. Such matrices appear in linear Robin
 * boundary operators, including thermal convection. Products `N_i N_j` have a
 * higher polynomial order than the interpolation used by ordinary surface
 * tractions, so this integration path selects a dedicated quadrature rule rather
 * than changing the established scheme returned by `Surface::integration_scheme`.
 *
 * The physical integration weight combines the reference quadrature weight with
 * the norm of the cross product of the two isoparametric surface tangents.
 *
 * @see Surface
 * @see SurfaceInterface::integrate_scalar_shape_matrix
 * @see ../../../../bc/robin/convection.h
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "surface.h"

namespace fem::model {

/**
 * Integrates a scalar-weighted consistent shape-function product matrix.
 *
 * The scalar callback is evaluated at every global quadrature-point position.
 * For each point, the routine evaluates the isoparametric shape vector, physical
 * surface Jacobian and complete physical area weight, and accumulates
 * `a(x) * N * N^T` into a fixed-size local matrix.
 *
 * Quadrature order is chosen from the surface interpolation order rather than
 * reusing ordinary load integration: linear three- and four-node surfaces use a
 * quadratic rule, while quadratic six- and eight-node surfaces use a quartic
 * rule. This is sufficient for the higher-order shape products while preserving
 * the existing pressure, traction and area integration behavior elsewhere.
 *
 * @param node_coords Global nodal coordinate field defining the physical surface.
 * @param field Scalar function evaluated in global coordinates.
 * @return Dense `N x N` local matrix in surface-connectivity ordering.
 */
template<Index N>
DynamicMatrix Surface<N>::integrate_scalar_shape_matrix(
    const Field&       node_coords,
    const ScalarField& field
) const {
    // Gather the physical nodal coordinates once because every quadrature point
    // uses the same isoparametric geometry.
    const auto coordinates = node_coords_global(node_coords);

    // Products N_i N_j need a higher-order rule than ordinary surface loads.
    // Keep this rule local to shape-matrix integration so mechanical pressure,
    // traction and area integration retain their established quadrature.
    const auto& scheme = []() -> const math::quadrature::Quadrature& {
        if constexpr (N == 3) {
            // Linear triangular interpolation produces quadratic shape products.
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_TRI,
                math::quadrature::ORDER_QUADRATIC
            };
            return q;
        } else if constexpr (N == 4) {
            // Bilinear quadrilateral shape products require quadratic integration
            // in each natural coordinate direction.
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_QUAD,
                math::quadrature::ORDER_QUADRATIC
            };
            return q;
        } else if constexpr (N == 6) {
            // Quadratic triangular interpolation generates quartic products.
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_TRI,
                math::quadrature::ORDER_QUARTIC
            };
            return q;
        } else {
            // Serendipity quadrilateral interpolation likewise requires the
            // dedicated quartic rule for its shape-function products.
            static const math::quadrature::Quadrature q{
                math::quadrature::DOMAIN_ISO_QUAD,
                math::quadrature::ORDER_QUARTIC
            };
            return q;
        }
    }();

    // Accumulate in the fixed compile-time surface size and return through the
    // dynamic interface expected by heterogeneous `SurfaceInterface` callers.
    StaticMatrix<N, N> matrix = StaticMatrix<N, N>::Zero();

    // Apply the selected quadrature rule over the complete natural element domain.
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        // Extract natural coordinates and the reference-domain quadrature weight.
        const auto point = scheme.get_point(ip);
        const Precision r = point.r;
        const Precision s = point.s;

        // Evaluate interpolation, physical tangents and global position at the
        // current natural quadrature point.
        const StaticMatrix<N, 1> shape = shape_function(r, s);
        const auto jac                   = jacobian(coordinates, r, s);
        const auto position              = interpolate(coordinates, r, s);

        // Convert the reference quadrature measure into physical surface area by
        // multiplying with the norm of the two tangent vectors' cross product.
        const Precision weighted_area = jac.col(0).cross(jac.col(1)).norm() * point.w;

        // Evaluate the supplied scalar coefficient in global coordinates.
        const Precision value = field(position);

        // Accumulate the consistent local matrix contribution a(x) N N^T dGamma.
        matrix.noalias() += value * weighted_area * (shape * shape.transpose());
    }

    return matrix;
}

} // namespace fem::model
