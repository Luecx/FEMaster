/**
 * @file surface3.h
 * @brief Declares the three-node linear triangular surface element.
 *
 * The element provides linear barycentric interpolation on the reference
 * triangle and delegates topology-independent geometry and integration to
 * `Surface<3>`.
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#pragma once

#include "surface.h"

#include <array>

namespace fem::model {

/**
 * @brief Three-node linear triangular finite-element surface.
 *
 * The element uses the natural triangular domain `r >= 0`, `s >= 0`, and
 * `r + s <= 1`. Its three nodes are located at `(0,0)`, `(1,0)` and `(0,1)`.
 * The boundary consists of three linear edges, and the element uses a
 * linear triangular quadrature rule.
 */
struct Surface3 : public Surface<3> {
    using Surface<3>::num_edges;
    using Surface<3>::num_nodes;
    using Surface<3>::num_nodes_per_edge;

    // Construction
    explicit Surface3(const std::array<ID, 3>& node_ids);

    // Linear interpolation on the reference triangle. The shape-function
    // derivatives define the constant physical tangent mapping of the element.
    StaticMatrix<3, 1> shape_function         (Precision r, Precision s) const override;
    StaticMatrix<3, 2> shape_derivative       (Precision r, Precision s) const override;
    StaticMatrix<3, 3> shape_second_derivative(Precision r, Precision s) const override;

    // Natural coordinates of the three corner nodes. The ordering must remain
    // consistent with the shape-function vector and the global connectivity.
    StaticMatrix<3, 2> node_coords_local() const override;

    // Projection onto the three linear boundary edges. The returned point is
    // expressed in the natural triangular coordinates used by this element.
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<3, 3>& node_coords) const override;

    // Check whether natural coordinates lie inside the closed reference
    // triangle, including points on its three boundary edges.
    bool in_bounds(const Vec2& local) const override;

    // Linear triangular quadrature rule used by all surface-field integration
    // routines inherited from the common surface implementation.
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model
