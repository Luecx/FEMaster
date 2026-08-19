/**
 * @file surface4.h
 * @brief Declares the four-node bilinear quadrilateral surface element.
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#pragma once

#include "surface.h"

namespace fem::model {

struct Surface4 : public Surface<4> {
    using Surface<4>::num_edges;
    using Surface<4>::num_nodes;
    using Surface<4>::num_nodes_per_edge;

    explicit Surface4(const std::array<ID, 4>& node_ids);

    Ptr copy() const override { return std::make_shared<Surface4>(*this); }

    StaticMatrix<4, 1> shape_function(Precision r, Precision s) const override;
    StaticMatrix<4, 2> shape_derivative(Precision r, Precision s) const override;
    StaticMatrix<4, 3> shape_second_derivative(Precision r, Precision s) const override;
    StaticMatrix<4, 2> node_coords_local() const override;
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<4, 3>& node_coords) const override;
    bool in_bounds(const Vec2& local) const override;
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model
