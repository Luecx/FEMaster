/**
 * @file surface6.h
 * @brief Declares the six-node quadratic triangular surface element.
 *
 * The element provides quadratic Lagrange interpolation on the reference
 * triangle and delegates topology-independent geometry and integration to
 * `Surface<6>`.
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
 * The element uses the natural triangular domain `r >= 0`, `s >= 0`, and
 * `r + s <= 1`. The first three nodes define the corners and the remaining
 * three nodes lie at the edge midpoints. The curved boundary is represented by
 * quadratic line elements, and the element uses a quadratic triangular
 * quadrature rule.
 */
struct Surface6 : public Surface<6> {
    using Surface<6>::num_edges;
    using Surface<6>::num_nodes;
    using Surface<6>::num_nodes_per_edge;

    // Construction
    explicit Surface6(const std::array<ID, 6>& node_ids);

    // Quadratic Lagrange interpolation on the reference triangle. The
    // derivatives describe the curved isoparametric geometry of the element.
    StaticMatrix<6, 1> shape_function         (Precision r, Precision s) const override;
    StaticMatrix<6, 2> shape_derivative       (Precision r, Precision s) const override;
    StaticMatrix<6, 3> shape_second_derivative(Precision r, Precision s) const override;

    // Natural coordinates of the three corner nodes followed by the three
    // midside nodes. This ordering is part of the element connectivity contract.
    StaticMatrix<6, 2> node_coords_local() const override;

    // Projection onto the three quadratic boundary edges. The returned point
    // is expressed in the natural triangular coordinates of this element.
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<6, 3>& node_coords) const override;

    // Check whether natural coordinates lie inside the closed reference
    // triangle, including points on its curved physical boundary.
    bool in_bounds(const Vec2& local) const override;

    // Quadratic triangular quadrature rule used by all surface-field
    // integration routines inherited from the common surface implementation.
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model
