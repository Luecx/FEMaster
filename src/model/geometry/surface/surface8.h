/**
* @file surface8.h
 * @brief Declares the eight-node quadratic quadrilateral surface element.
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
 * The element uses the natural quadrilateral domain
 * `-1 <= r <= 1` and `-1 <= s <= 1`.
 */
struct Surface8 : public Surface<8> {
    using Surface<8>::num_edges;
    using Surface<8>::num_nodes;
    using Surface<8>::num_nodes_per_edge;

    // Construction
    explicit Surface8(const std::array<ID, 8>& node_ids);

    // Shape functions
    StaticMatrix<8, 1> shape_function         (Precision r, Precision s) const override;
    StaticMatrix<8, 2> shape_derivative       (Precision r, Precision s) const override;
    StaticMatrix<8, 3> shape_second_derivative(Precision r, Precision s) const override;

    // local values of node coordinates in natural coordinates
    StaticMatrix<8, 2> node_coords_local() const override;

    // Boundary projection and local domain
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<8, 3>& node_coords) const override;

    // check if a local value (r,s) is inside or not
    bool in_bounds(const Vec2& local) const override;

    // Numerical integration
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model