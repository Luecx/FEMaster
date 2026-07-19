/**
* @file surface6.h
 * @brief Declares the six-node quadratic triangular surface element.
 *
 * @author Finn Eggers
 * @date 01.10.2024
 */

#pragma once

#include "surface.h"

#include <array>

namespace fem::model {

/**
 * @brief Six-node quadratic triangular finite-element surface.
 *
 * The element uses the natural triangular domain
 * `r >= 0`, `s >= 0`, and `r + s <= 1`.
 */
struct Surface6 : public Surface<6> {
    using Surface<6>::num_edges;
    using Surface<6>::num_nodes;
    using Surface<6>::num_nodes_per_edge;

    // Construction
    explicit Surface6(const std::array<ID, 6>& node_ids);

    // Shape functions
    StaticMatrix<6, 1> shape_function         (Precision r, Precision s) const override;
    StaticMatrix<6, 2> shape_derivative       (Precision r, Precision s) const override;
    StaticMatrix<6, 3> shape_second_derivative(Precision r, Precision s) const override;

    // local values of node coordinates in natural coordinates
    StaticMatrix<6, 2> node_coords_local() const override;

    // Boundary projection and local domain
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<6, 3>& node_coords) const override;

    // check if a local value (r,s) is inside or not
    bool in_bounds(const Vec2& local) const override;

    // Numerical integration
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model