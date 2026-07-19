/**
 * @file surface8.h
 * @brief Declares the eight-node quadratic quadrilateral surface element.
 *
 * The element provides serendipity interpolation on the reference square and
 * delegates topology-independent geometry and integration to `Surface<8>`.
 *
 * @author Finn Eggers
 * @date 01.10.2024
 */

#pragma once

#include "surface.h"

#include <array>

namespace fem::model {

/**
 * @brief Eight-node serendipity quadrilateral finite-element surface.
 *
 * The element uses the natural quadrilateral domain `-1 <= r <= 1` and
 * `-1 <= s <= 1`. The first four nodes define the corners and the remaining
 * four nodes lie at the edge midpoints. The curved boundary is represented by
 * quadratic line elements, and the element uses a quadratic quadrilateral
 * quadrature rule.
 */
struct Surface8 : public Surface<8> {
    using Surface<8>::num_edges;
    using Surface<8>::num_nodes;
    using Surface<8>::num_nodes_per_edge;

    // Construction
    explicit Surface8(const std::array<ID, 8>& node_ids);

    // Serendipity interpolation on the reference square. The first and second
    // derivatives describe the curved isoparametric geometry of the element.
    StaticMatrix<8, 1> shape_function         (Precision r, Precision s) const override;
    StaticMatrix<8, 2> shape_derivative       (Precision r, Precision s) const override;
    StaticMatrix<8, 3> shape_second_derivative(Precision r, Precision s) const override;

    // Natural coordinates of the four corner nodes followed by the four
    // midside nodes. This ordering is part of the element connectivity contract.
    StaticMatrix<8, 2> node_coords_local() const override;

    // Projection onto the four quadratic boundary edges. The returned point
    // is expressed in the natural quadrilateral coordinates of this element.
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<8, 3>& node_coords) const override;

    // Check whether natural coordinates lie inside the closed reference
    // square, including points on its curved physical boundary.
    bool in_bounds(const Vec2& local) const override;

    // Quadratic quadrilateral quadrature rule used by all surface-field
    // integration routines inherited from the common surface implementation.
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model
