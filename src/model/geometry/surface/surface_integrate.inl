/**
 * @file surface_integrate.h
 * @brief Implements field integration for finite-element surfaces.
 */

#pragma once

#include "surface.h"

#include <cmath>
#include <functional>

namespace fem::model {

template<Index N>
Precision Surface<N>::area(const Field& node_coords) const {
    // the surface area is obtained by integrating the constant scalar field one
    // over the complete physical surface
    return integrate_scalar_field(node_coords, [](const Vec3&) {
        return Precision(1);
    });
}

template<Index N>
Precision Surface<N>::integrate_scalar_field(
    const Field&       node_coords,
    const ScalarField& field
) const {
    // collect the global coordinates of all nodes belonging to this surface
    const auto coordinates = node_coords_global(node_coords);

    // the quadrature scheme integrates in the natural coordinate system of the
    // element. therefore the field value has to be multiplied by the physical
    // surface jacobian at every integration point.
    std::function<Precision(Precision, Precision, Precision)> integrand =
        [&](Precision r, Precision s, Precision) {
            const auto jac      = jacobian(coordinates, r, s);
            const auto position = interpolate(coordinates, r, s);

            // the cross product of both surface tangent vectors transforms an
            // infinitesimal natural-coordinate area into physical surface area
            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            return surface_jacobian * field(position);
        };

    // apply the integration points and weights of the complete surface element
    return integration_scheme().integrate(integrand);
}

template<Index N>
Vec3 Surface<N>::integrate_vector_field(
    const Field&    node_coords,
    const VecField& field
) const {
    // collect the global coordinates of all nodes belonging to this surface
    const auto coordinates = node_coords_global(node_coords);

    // vector-valued integration uses the same geometric area transformation as
    // scalar integration; only the returned field value differs
    std::function<Vec3(Precision, Precision, Precision)> integrand =
        [&](Precision r, Precision s, Precision) {
            const auto jac      = jacobian(coordinates, r, s);
            const auto position = interpolate(coordinates, r, s);

            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            return surface_jacobian * field(position);
        };

    return integration_scheme().integrate(integrand);
}

template<Index N>
void Surface<N>::integrate_vector_field(
    const Field&    node_coords,
    Field&          target,
    const VecField& field
) const {
    // collect the global coordinates once because they are reused at every
    // integration point
    const auto coordinates = node_coords_global(node_coords);
    const auto& scheme      = integration_scheme();

    // this overload does not return the total vector integral. instead it
    // distributes the integrated field consistently onto the surface nodes
    // using the shape functions.
    for (Index local_ip = 0; local_ip < scheme.count(); ++local_ip) {
        const auto point = scheme.get_point(local_ip);

        const Precision r = point.r;
        const Precision s = point.s;
        const Precision w = point.w;

        // the shape functions distribute the contribution at the current
        // integration point onto the individual surface nodes
        const StaticMatrix<N, 1> shape = shape_function(r, s);

        // evaluate the physical position and the local surface area scaling
        const auto jac      = jacobian(coordinates, r, s);
        const auto position = interpolate(coordinates, r, s);

        // include the quadrature weight directly in the differential area
        const Precision weighted_area = jac.col(0).cross(jac.col(1)).norm() * w;

        // evaluate the prescribed vector field in physical coordinates
        const Vec3 value = field(position);

        // assemble the consistent nodal contribution
        // f_i += N_i(r,s) * field(x(r,s)) * dA
        for (Index local_id = 0; local_id < N; ++local_id) {
            for (Dim dof = 0; dof < 3; ++dof) {
                target(nodeIds[local_id], dof) +=
                    shape[local_id] * value[dof] * weighted_area;
            }
        }
    }
}

template<Index N>
Mat3 Surface<N>::integrate_tensor_field(
    const Field&    node_coords,
    const TenField& field
) const {
    // collect the global coordinates of all nodes belonging to this surface
    const auto coordinates = node_coords_global(node_coords);

    // tensor-valued integration follows the same local-to-global area
    // transformation as scalar and vector integration
    std::function<Mat3(Precision, Precision, Precision)> integrand =
        [&](Precision r, Precision s, Precision) {
            const auto jac      = jacobian(coordinates, r, s);
            const auto position = interpolate(coordinates, r, s);

            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            return surface_jacobian * field(position);
        };

    return integration_scheme().integrate(integrand);
}

template<Index N>
void Surface<N>::integrate_triangular(
    const Field&                                                    node_coords,
    const Polygon&                                                  polygon,
    const math::quadrature::Quadrature&                             scheme,
    const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand
) const {
    // the supplied polygon may extend beyond the natural coordinate domain of
    // the surface. clip it first so that only the valid overlap is integrated.
    const auto intersection = polygon.intersection(this->local_domain_polygon());

    // fewer than three points cannot form a finite integration region
    if (intersection.size() < 3) {
        return;
    }

    // an intersection may contain several points while still being numerically
    // degenerate due to coincident or nearly collinear vertices
    if (intersection.area() <= Precision(1e-12)) {
        return;
    }

    // collect the global nodal coordinates once for all generated triangles
    // and quadrature points
    const auto coordinates = node_coords_global(node_coords);

    // the convex intersection polygon is divided into triangles that all share
    // the first polygon point
    const Vec2 origin = intersection[0];

    // triangulate the intersection polygon using the triangle fan
    // (p0,p1,p2), (p0,p2,p3), ..., (p0,p[n-2],p[n-1])
    for (std::size_t i = 1; i + 1 < intersection.size(); ++i) {
        const Vec2 edge_r = intersection[i]     - origin;
        const Vec2 edge_s = intersection[i + 1] - origin;

        // the determinant is the area scaling from the reference triangle onto
        // the current triangle in the natural coordinates of the surface
        const Precision triangle_jacobian =
            std::abs(edge_r(0) * edge_s(1) - edge_r(1) * edge_s(0));

        // duplicate or collinear fan points create a zero-area triangle and
        // therefore do not contribute to the integral
        if (triangle_jacobian <= Precision(1e-12)) {
            continue;
        }

        // apply the supplied triangular quadrature rule to the current triangle
        for (Index q = 0; q < scheme.count(); ++q) {
            const auto point = scheme.get_point(q);

            const Precision r = point.r;
            const Precision s = point.s;
            const Precision w = point.w;

            // map the quadrature point from the reference triangle onto the
            // current triangle in the natural coordinate system of the surface
            const Vec2 local = origin + r * edge_r + s * edge_s;

            // map the natural surface coordinates further into physical space
            const Vec3 global = interpolate(coordinates, local(0), local(1));

            // the finite-element surface may be curved, so its physical area
            // scaling has to be evaluated at every quadrature point
            const auto jac = jacobian(coordinates, local(0), local(1));
            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            // combine the quadrature weight, the reference-to-local triangle
            // mapping and the local-to-global surface mapping
            const Precision weight = w * triangle_jacobian * surface_jacobian;

            // provide the local position, physical position and complete
            // physical integration weight to the caller
            integrand(local, global, weight);
        }
    }
}

} // namespace fem::model
