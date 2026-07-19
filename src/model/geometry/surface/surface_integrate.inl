/**
 * @file surface_integrate.inl
 * @brief Implements numerical integration of fields over finite-element surfaces.
 *
 * Scalar, vector and tensor fields are evaluated at quadrature points in
 * natural coordinates and weighted with the physical surface Jacobian. The
 * polygon-restricted routine additionally clips the natural integration domain
 * and integrates the overlap by a triangle fan.
 *
 * @see Surface
 * @see SurfacePolygon
 */

#pragma once

#include "surface.h"

#include <cmath>
#include <functional>

namespace fem::model {

/**
 * Computes the physical area of the surface.
 *
 * The area is obtained by integrating the constant scalar field one over the
 * complete natural element domain. The integration routine applies the
 * physical surface Jacobian at each quadrature point.
 *
 * @param node_coords Global nodal coordinate field.
 *
 * @return Physical surface area.
 */
template<Index N>
Precision Surface<N>::area(const Field& node_coords) const {
    // Integrate the constant unit field over the complete physical surface
    return integrate_scalar_field(node_coords, [](const Vec3&) {
        return Precision(1);
    });
}

/**
 * Integrates a scalar field over the complete physical surface.
 *
 * The quadrature rule is defined in natural coordinates. At every quadrature
 * point, the scalar field is evaluated at the interpolated global position and
 * multiplied by the norm of the cross product of the two physical tangent
 * vectors. This converts the natural-coordinate measure into physical area.
 *
 * @param node_coords Global nodal coordinate field.
 * @param field Scalar field evaluated at global positions.
 *
 * @return Integral of the scalar field over the physical surface.
 */
template<Index N>
Precision Surface<N>::integrate_scalar_field(
    const Field&       node_coords,
    const ScalarField& field
) const {
    // Gather the global coordinates once because they are reused at every
    // quadrature point
    const auto coordinates = node_coords_global(node_coords);

    // Convert the natural-coordinate quadrature measure to physical surface
    // area at every integration point
    std::function<Precision(Precision, Precision, Precision)> integrand =
        [&](Precision r, Precision s, Precision) {
            const auto jac      = jacobian(coordinates, r, s);
            const auto position = interpolate(coordinates, r, s);

            // The cross product of the surface tangents supplies the physical
            // area scaling of the isoparametric mapping
            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            return surface_jacobian * field(position);
        };

    // Apply the quadrature rule over the complete surface element
    return integration_scheme().integrate(integrand);
}

/**
 * Integrates a vector field over the complete physical surface.
 *
 * The vector field is evaluated at the interpolated global position of each
 * quadrature point. Its value is multiplied by the physical surface Jacobian
 * before the element quadrature rule accumulates the result.
 *
 * @param node_coords Global nodal coordinate field.
 * @param field Vector field evaluated at global positions.
 *
 * @return Integral of the vector field over the physical surface.
 */
template<Index N>
Vec3 Surface<N>::integrate_vector_field(
    const Field&    node_coords,
    const VecField& field
) const {
    // Gather the global coordinates once because they are reused at every
    // quadrature point
    const auto coordinates = node_coords_global(node_coords);

    // Use the same physical area transformation as scalar integration while
    // retaining the vector-valued field result
    std::function<Vec3(Precision, Precision, Precision)> integrand =
        [&](Precision r, Precision s, Precision) {
            const auto jac      = jacobian(coordinates, r, s);
            const auto position = interpolate(coordinates, r, s);

            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            return surface_jacobian * field(position);
        };

    return integration_scheme().integrate(integrand);
}

/**
 * Integrates a vector field and assembles its consistent nodal contribution.
 *
 * At each quadrature point, the vector field is evaluated in global
 * coordinates and distributed to the surface nodes with the shape functions.
 * The contribution includes both the physical surface Jacobian and the
 * quadrature weight, so the target field receives the complete physical
 * integral.
 *
 * @param node_coords Global nodal coordinate field.
 * @param target Mutable global field receiving the nodal contributions.
 * @param field Vector field evaluated at global positions.
 */
template<Index N>
void Surface<N>::integrate_vector_field(
    const Field&    node_coords,
    Field&          target,
    const VecField& field
) const {
    // Gather the global coordinates once because they are reused at every
    // integration point
    const auto coordinates = node_coords_global(node_coords);
    const auto& scheme      = integration_scheme();

    // Distribute the integrated field consistently onto the surface nodes
    // instead of returning one total vector integral
    for (Index local_ip = 0; local_ip < scheme.count(); ++local_ip) {
        const auto point = scheme.get_point(local_ip);

        const Precision r = point.r;
        const Precision s = point.s;
        const Precision w = point.w;

        // Evaluate the shape functions used to distribute the current
        // integration-point contribution onto the surface nodes
        const StaticMatrix<N, 1> shape = shape_function(r, s);

        // Evaluate the physical position and local surface area scaling
        const auto jac      = jacobian(coordinates, r, s);
        const auto position = interpolate(coordinates, r, s);

        // Include the quadrature weight in the differential physical area
        const Precision weighted_area = jac.col(0).cross(jac.col(1)).norm() * w;

        // Evaluate the prescribed vector field in physical coordinates
        const Vec3 value = field(position);

        // Assemble the consistent nodal contribution
        // f_i += N_i(r,s) * field(x(r,s)) * dA
        for (Index local_id = 0; local_id < N; ++local_id) {
            for (Dim dof = 0; dof < 3; ++dof) {
                target(nodeIds[local_id], dof) +=
                    shape[local_id] * value[dof] * weighted_area;
            }
        }
    }
}

/**
 * Integrates a tensor field over the complete physical surface.
 *
 * Tensor values are evaluated at interpolated global positions and weighted
 * with the physical surface Jacobian before applying the element quadrature
 * rule.
 *
 * @param node_coords Global nodal coordinate field.
 * @param field Tensor field evaluated at global positions.
 *
 * @return Integral of the tensor field over the physical surface.
 */
template<Index N>
Mat3 Surface<N>::integrate_tensor_field(
    const Field&    node_coords,
    const TenField& field
) const {
    // Gather the global coordinates once because they are reused at every
    // quadrature point
    const auto coordinates = node_coords_global(node_coords);

    // Apply the same natural-to-physical area transformation as scalar and
    // vector integration
    std::function<Mat3(Precision, Precision, Precision)> integrand =
        [&](Precision r, Precision s, Precision) {
            const auto jac      = jacobian(coordinates, r, s);
            const auto position = interpolate(coordinates, r, s);

            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            return surface_jacobian * field(position);
        };

    return integration_scheme().integrate(integrand);
}

/**
 * Integrates the overlap of a polygon with the natural surface domain.
 *
 * The supplied polygon is clipped against the natural element domain. The
 * valid overlap is decomposed into a triangle fan, and the supplied
 * triangular quadrature rule is applied to every fan triangle. Each callback
 * receives the natural coordinates, global coordinates and complete physical
 * integration weight. That weight combines the quadrature weight, the
 * natural-triangle mapping and the physical surface Jacobian.
 *
 * Degenerate intersections and zero-area fan triangles are ignored.
 *
 * @param node_coords Global nodal coordinate field.
 * @param polygon Integration polygon in natural surface coordinates.
 * @param scheme Triangular quadrature rule on the reference triangle.
 * @param integrand Callback receiving natural coordinates, global coordinates
 *                  and the complete physical integration weight.
 */
template<Index N>
void Surface<N>::integrate_triangular(
    const Field&                                                    node_coords,
    const Polygon&                                                  polygon,
    const math::quadrature::Quadrature&                             scheme,
    const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand
) const {
    // Clip the supplied polygon against the valid natural-coordinate domain
    const auto intersection = polygon.intersection(this->local_domain_polygon());

    // Reject polygons that cannot form a finite integration region
    if (intersection.size() < 3) {
        return;
    }

    // Reject numerically degenerate intersections with coincident or nearly
    // collinear vertices
    if (intersection.area() <= Precision(1e-12)) {
        return;
    }

    // Gather the global nodal coordinates once for all triangles and
    // quadrature points
    const auto coordinates = node_coords_global(node_coords);

    // Decompose the convex intersection polygon into triangles sharing its
    // first point
    const Vec2 origin = intersection[0];

    // Apply the triangle-fan decomposition
    for (std::size_t i = 1; i + 1 < intersection.size(); ++i) {
        const Vec2 edge_r = intersection[i]     - origin;
        const Vec2 edge_s = intersection[i + 1] - origin;

        // Compute the area scaling from the reference triangle to the current
        // triangle in natural coordinates
        const Precision triangle_jacobian =
            std::abs(edge_r(0) * edge_s(1) - edge_r(1) * edge_s(0));

        // Ignore zero-area fan triangles caused by duplicate or collinear
        // points
        if (triangle_jacobian <= Precision(1e-12)) {
            continue;
        }

        // Apply the supplied triangular quadrature rule to the current triangle
        for (Index q = 0; q < scheme.count(); ++q) {
            const auto point = scheme.get_point(q);

            const Precision r = point.r;
            const Precision s = point.s;
            const Precision w = point.w;

            // Map the quadrature point from the reference triangle into the
            // current triangle in natural coordinates
            const Vec2 local = origin + r * edge_r + s * edge_s;

            // Map the natural surface coordinates into physical space
            const Vec3 global = interpolate(coordinates, local(0), local(1));

            // Evaluate the physical area scaling at the current point because
            // the finite-element surface may be curved
            const auto jac = jacobian(coordinates, local(0), local(1));
            const Precision surface_jacobian = jac.col(0).cross(jac.col(1)).norm();

            // Combine the quadrature weight, the natural triangle mapping and
            // the physical surface mapping
            const Precision weight = w * triangle_jacobian * surface_jacobian;

            // Pass the natural position, physical position and complete
            // physical integration weight to the caller
            integrand(local, global, weight);
        }
    }
}

} // namespace fem::model
