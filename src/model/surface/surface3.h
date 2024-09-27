#pragma once

#include "surface.h"

#include <array>
#include <tuple>

namespace fem::model {

/******************************************************************************
 * @class Surface3
 * @brief Derived class for a triangular surface element with 3 nodes (S3D3).
 ******************************************************************************/
struct Surface3 : public SurfaceInterface<3> {

    Surface3(const std::array<ID, 3>& pNodeIds)
        : SurfaceInterface<3>(pNodeIds) {}

    StaticMatrix<3, 1> shape_function(Precision r, Precision s) const override {
        StaticMatrix<3, 1> N;
        N(0, 0) = 1 - r - s;  // Shape function for node 1
        N(1, 0) = r;          // Shape function for node 2
        N(2, 0) = s;          // Shape function for node 3
        return N;
    }

    StaticMatrix<3, 2> shape_derivative(Precision r, Precision s) const override {
        StaticMatrix<3, 2> dN;
        dN(0, 0) = -1;  dN(0, 1) = -1;  // dN1/dr, dN1/ds
        dN(1, 0) = 1;   dN(1, 1) = 0;   // dN2/dr, dN2/ds
        dN(2, 0) = 0;   dN(2, 1) = 1;   // dN3/dr, dN3/ds
        return dN;
    }

    StaticMatrix<3, 2> node_coords_local() const override {
        StaticMatrix<3, 2> local_coords;
        local_coords << 0.0, 0.0,   // Node 1
                        1.0, 0.0,   // Node 2
                        0.0, 1.0;   // Node 3
        return local_coords;
    }

    StaticVector<2> global_to_local(const StaticVector<3>& global, const NodeData& node_coords_system, bool clip) const override {
        StaticVector<2> local;
        auto node_coords = this->node_coords_global(node_coords_system);

        StaticVector<3> v0 = node_coords.row(2) - node_coords.row(0);  // Node 3 - Node 1
        StaticVector<3> v1 = node_coords.row(1) - node_coords.row(0);  // Node 2 - Node 1
        StaticVector<3> v2 = global - node_coords.row(0).transpose();  // From Node 1 to point

        Precision d00 = v0.dot(v0);
        Precision d01 = v0.dot(v1);
        Precision d11 = v1.dot(v1);
        Precision d20 = v2.dot(v0);
        Precision d21 = v2.dot(v1);
        Precision denom = d00 * d11 - d01 * d01;

        local(0) = (d00 * d21 - d01 * d20) / denom;  // r
        local(1) = (d11 * d20 - d01 * d21) / denom;  // s

        if (clip) {
            if (local(0) < 0 || local(1) < 0 || local(0) + local(1) > 1) {
                local = closest_point_on_triangle_edges(global, node_coords);
            }
        }

        return local;
    }

    StaticVector<2> closest_point_on_triangle_edges(const StaticVector<3>& global, const StaticMatrix<3, 3>& node_coords) const {
        auto res1 = closest_point_on_edge(global, node_coords.row(0), node_coords.row(1));
        auto res2 = closest_point_on_edge(global, node_coords.row(1), node_coords.row(2));
        auto res3 = closest_point_on_edge(global, node_coords.row(2), node_coords.row(0));

        Precision d1 = (std::get<0>(res1) - global).squaredNorm();
        Precision d2 = (std::get<0>(res2) - global).squaredNorm();
        Precision d3 = (std::get<0>(res3) - global).squaredNorm();

        Precision best_dist = std::min({d1, d2, d3});

        if (d1 == best_dist) return {std::get<1>(res1), 0};
        if (d2 == best_dist) return {1 - std::get<1>(res2), std::get<1>(res2)};
        return {0, 1-std::get<1>(res3)};
    }

    std::tuple<StaticVector<3>, Precision> closest_point_on_edge(
        const StaticVector<3>& global, const StaticVector<3>& p0, const StaticVector<3>& p1) const {
        StaticVector<3> edge = p1 - p0;
        StaticVector<3> v2 = global - p0;

        Precision t = v2.dot(edge) / edge.squaredNorm();

        t = std::max(0.0, std::min(1.0, t));
        StaticVector<3> closest_point = p0 + t * edge;
        return std::make_tuple(closest_point, t);
    }

    const quadrature::Quadrature& integration_scheme() const override {
        static const quadrature::Quadrature quad{quadrature::DOMAIN_ISO_TRI, quadrature::ORDER_LINEAR};
        return quad;
    }
};

}  // namespace fem::model
