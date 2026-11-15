/**
 * @file surface_load.inl
 * @brief Implements consistent surface loads and their analytic load tangents.
 *
 * Pressure and traction loads are integrated over the surface configuration
 * supplied by the caller. Dead loads pass the reference positions and request
 * no tangent. Follower loads pass the current positions and additionally
 * receive the exact derivative of the quadrature-based nodal force with
 * respect to the nodal positions.
 *
 * With the covariant tangent vectors a_r = x_,r and a_s = x_,s, the
 * unnormalized normal c = a_r x a_s carries both the surface orientation and
 * the area measure, dGamma = |c| dr ds. Its derivative with respect to the
 * position of node k is
 *
 *     dc/dx_k = N_k,s [a_r]x - N_k,r [a_s]x,
 *
 * which is the common building block of all surface load tangents.
 *
 * Results use surface-connectivity ordering with three translational
 * components per node. Scattering into global fields and equation indices
 * belongs to the calling load.
 *
 * @see Surface
 * @see bc::PLoad
 * @see bc::DLoad
 *
 * @author Finn Eggers
 * @date 15.11.2026
 */

#pragma once

#include "surface.h"

#include "../../../core/logging.h"
#include "../../../math/vec_util.h"

#include <cmath>

namespace fem::model {

/**
 * Integrates a pressure acting against the geometric surface normal.
 *
 * The consistent nodal force of a pressure p acting against the orientation of
 * the surface is
 *
 *     f_i = -p integral N_i n dGamma = -p sum_q w_q N_i(q) c(q),
 *
 * with the unnormalized normal c = a_r x a_s. No normalization is required
 * because |c| is exactly the area measure.
 *
 * If requested, the tangent receives the analytic derivative with respect to
 * the nodal positions of the supplied configuration,
 *
 *     df_i/dx_k = -p sum_q w_q N_i(q) (N_k,s [a_r]x - N_k,r [a_s]x).
 *
 * This derivative is generally nonsymmetric. It becomes symmetric only for
 * closed surfaces or suitable boundary conditions.
 *
 * @param node_coords Nodal positions defining the loaded configuration.
 * @param pressure Effective pressure including the caller's multiplier.
 * @param force Local nodal force vector of size 3N, overwritten.
 * @param tangent Optional local derivative df/dx of size 3N x 3N, overwritten.
 */
template<Index N>
void Surface<N>::integrate_pressure(
    const Field&   node_coords,
    Precision      pressure,
    DynamicVector& force,
    DynamicMatrix* tangent
) const {
    // Gather the nodal positions of the loaded configuration once
    const StaticMatrix<N, 3> coordinates = node_coords_global(node_coords);
    const auto&              scheme      = integration_scheme();

    StaticVector<3 * N>        f = StaticVector<3 * N>::Zero();
    StaticMatrix<3 * N, 3 * N> k = StaticMatrix<3 * N, 3 * N>::Zero();

    // Apply the surface quadrature rule used by all surface loads so the
    // tangent is the exact derivative of the discrete residual
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        // Evaluate the interpolation and its natural derivatives
        const StaticMatrix<N, 1> shape = shape_function  (point.r, point.s);
        const StaticMatrix<N, 2> deriv = shape_derivative(point.r, point.s);

        // Covariant tangent vectors and unnormalized normal of the configuration
        const StaticMatrix<3, 2> jac = coordinates.transpose() * deriv;
        const Vec3               a_r = jac.col(0);
        const Vec3               a_s = jac.col(1);
        const Vec3               c   = a_r.cross(a_s);

        // Pressure traction including the area measure and quadrature weight
        const Vec3 traction = -pressure * point.w * c;

        for (Index i = 0; i < N; ++i) {
            f.template segment<3>(3 * i) += shape(i) * traction;
        }

        if (tangent == nullptr) continue;

        // Linearize the unnormalized normal with respect to every surface node
        const Mat3 skew_r = math::skew(a_r);
        const Mat3 skew_s = math::skew(a_s);

        for (Index node_k = 0; node_k < N; ++node_k) {
            const Mat3 dc_dx = deriv(node_k, 1) * skew_r - deriv(node_k, 0) * skew_s;
            const Mat3 dt_dx = -pressure * point.w * dc_dx;

            for (Index i = 0; i < N; ++i) {
                k.template block<3, 3>(3 * i, 3 * node_k) += shape(i) * dt_dx;
            }
        }
    }

    // Return the local results through the dynamic interface
    force = f;
    if (tangent != nullptr) *tangent = k;
}

/**
 * Integrates a traction defined per unit area of the supplied configuration.
 *
 * The consistent nodal force of the traction field t(x) is
 *
 *     f_i = integral N_i t(x) dGamma = sum_q w_q N_i(q) t(x_q) |c(q)|.
 *
 * The traction direction is fixed in space or given by a position-dependent
 * basis. If requested, the analytic derivative with respect to the nodal
 * positions combines the change of the area measure and of the traction:
 *
 *     d|c|/dx_k   = n^T (N_k,s [a_r]x - N_k,r [a_s]x),
 *     dt/dx_k     = (dt/dx)(x_q) N_k(q),
 *     df_i/dx_k   = sum_q w_q N_i(q) (t n^T dc/dx_k + |c| dt/dx N_k).
 *
 * A missing traction gradient denotes a spatially constant traction. A
 * degenerate surface without a finite normal cannot carry a traction and is
 * reported as an error.
 *
 * @param node_coords Nodal positions defining the loaded configuration.
 * @param traction Global traction per unit area including the caller's multiplier.
 * @param traction_gradient Optional spatial gradient dt/dx of the traction.
 * @param force Local nodal force vector of size 3N, overwritten.
 * @param tangent Optional local derivative df/dx of size 3N x 3N, overwritten.
 */
template<Index N>
void Surface<N>::integrate_traction(
    const Field&    node_coords,
    const VecField& traction,
    const TenField* traction_gradient,
    DynamicVector&  force,
    DynamicMatrix*  tangent
) const {
    // Gather the nodal positions of the loaded configuration once
    const StaticMatrix<N, 3> coordinates = node_coords_global(node_coords);
    const auto&              scheme      = integration_scheme();

    StaticVector<3 * N>        f = StaticVector<3 * N>::Zero();
    StaticMatrix<3 * N, 3 * N> k = StaticMatrix<3 * N, 3 * N>::Zero();

    // Apply the surface quadrature rule used by all surface loads so the
    // tangent is the exact derivative of the discrete residual
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        // Evaluate the interpolation and its natural derivatives
        const StaticMatrix<N, 1> shape = shape_function  (point.r, point.s);
        const StaticMatrix<N, 2> deriv = shape_derivative(point.r, point.s);

        // Covariant tangent vectors, area measure and physical position
        const StaticMatrix<3, 2> jac      = coordinates.transpose() * deriv;
        const Vec3               a_r      = jac.col(0);
        const Vec3               a_s      = jac.col(1);
        const Vec3               c        = a_r.cross(a_s);
        const Precision          area     = c.norm();
        const Vec3               position = coordinates.transpose() * shape;

        logging::error(area > Precision(0) && std::isfinite(area),
            "SURFACE: traction integration requires a non-degenerate surface");

        // Traction including the area measure and quadrature weight
        const Vec3 t = traction(position);

        for (Index i = 0; i < N; ++i) {
            f.template segment<3>(3 * i) += shape(i) * area * point.w * t;
        }

        if (tangent == nullptr) continue;

        // Linearize the area measure and the position-dependent traction
        const Vec3 normal = c / area;
        const Mat3 skew_r = math::skew(a_r);
        const Mat3 skew_s = math::skew(a_s);
        const Mat3 grad_t = traction_gradient ? (*traction_gradient)(position) : Mat3(Mat3::Zero());

        for (Index node_k = 0; node_k < N; ++node_k) {
            const Mat3 dc_dx   = deriv(node_k, 1) * skew_r - deriv(node_k, 0) * skew_s;
            const Mat3 d_area  = t * (normal.transpose() * dc_dx);
            const Mat3 d_field = area * shape(node_k) * grad_t;
            const Mat3 d_load  = point.w * (d_area + d_field);

            for (Index i = 0; i < N; ++i) {
                k.template block<3, 3>(3 * i, 3 * node_k) += shape(i) * d_load;
            }
        }
    }

    // Return the local results through the dynamic interface
    force = f;
    if (tangent != nullptr) *tangent = k;
}

} // namespace fem::model
