/**
 * @file surface4.h
 * @brief Declares the four-node bilinear quadrilateral surface element.
 *
 * The element provides bilinear corner interpolation on the reference square
 * and delegates topology-independent geometry and integration to `Surface<4>`.
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#pragma once

#include "surface.h"

#include <array>

namespace fem::model {

/**
 * @brief Four-node bilinear quadrilateral finite-element surface.
 *
 * The element uses the natural quadrilateral domain `-1 <= r <= 1` and
 * `-1 <= s <= 1`. Its four nodes are ordered counter-clockwise around the
 * reference square. The boundary consists of four linear edges, and the
 * element uses a linear quadrilateral quadrature rule.
 */
struct Surface4 : public Surface<4> {
    using Surface<4>::num_edges;
    using Surface<4>::num_nodes;
    using Surface<4>::num_nodes_per_edge;

    // Construction
    explicit Surface4(const std::array<ID, 4>& node_ids);

    // Bilinear interpolation on the reference square. The first derivatives
    // provide the varying physical tangent mapping of the element.
    StaticMatrix<4, 1> shape_function         (Precision r, Precision s) const override;
    StaticMatrix<4, 2> shape_derivative       (Precision r, Precision s) const override;
    StaticMatrix<4, 3> shape_second_derivative(Precision r, Precision s) const override;

    // Natural coordinates of the four corner nodes. The ordering follows the
    // counter-clockwise boundary orientation of the reference square.
    StaticMatrix<4, 2> node_coords_local() const override;

    // Projection onto the four linear boundary edges. The returned point is
    // expressed in the natural quadrilateral coordinates of this element.
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<4, 3>& node_coords) const override;

    // Check whether natural coordinates lie inside the closed reference
    // square, including points on its boundary.
    bool in_bounds(const Vec2& local) const override;

    // Linear quadrilateral quadrature rule used by all surface-field
    // integration routines inherited from the common surface implementation.
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model
