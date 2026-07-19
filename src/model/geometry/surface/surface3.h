/**
* @file surface3.h
 * @brief Declares the three-node linear triangular surface element.
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
 * The element uses the natural triangular domain
 * `r >= 0`, `s >= 0`, and `r + s <= 1`.
 */
struct Surface3 : public Surface<3> {
    using Surface<3>::num_edges;
    using Surface<3>::num_nodes;
    using Surface<3>::num_nodes_per_edge;

    // Construction
    explicit Surface3(const std::array<ID, 3>& node_ids);

    // Shape functions
    StaticMatrix<3, 1> shape_function         (Precision r, Precision s) const override;
    StaticMatrix<3, 2> shape_derivative       (Precision r, Precision s) const override;
    StaticMatrix<3, 3> shape_second_derivative(Precision r, Precision s) const override;

    // local values of node coordinates in natural coordinates
    StaticMatrix<3, 2> node_coords_local() const override;

    // Boundary projection and local domain
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<3, 3>& node_coords) const override;

    // check if a local value (r,s) is inside or not.
    bool in_bounds(const Vec2& local) const override;

    // Numerical integration
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model