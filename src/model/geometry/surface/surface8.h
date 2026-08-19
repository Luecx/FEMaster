/**
 * @file surface8.h
 * @brief Declares the eight-node quadratic quadrilateral surface element.
 *
 * @author Finn Eggers
 * @date 01.10.2024
 */

#pragma once

#include "surface.h"

namespace fem::model {

struct Surface8 : public Surface<8> {
    using Surface<8>::num_edges;
    using Surface<8>::num_nodes;
    using Surface<8>::num_nodes_per_edge;

    explicit Surface8(const std::array<ID, 8>& node_ids);

    Ptr copy() const override { return std::make_shared<Surface8>(*this); }

    StaticMatrix<8, 1> shape_function(Precision r, Precision s) const override;
    StaticMatrix<8, 2> shape_derivative(Precision r, Precision s) const override;
    StaticMatrix<8, 3> shape_second_derivative(Precision r, Precision s) const override;
    StaticMatrix<8, 2> node_coords_local() const override;
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<8, 3>& node_coords) const override;
    bool in_bounds(const Vec2& local) const override;
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace fem::model
