#pragma once

#include "surface.h"
#include "../line/line2b.h"

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

    StaticMatrix<3, 3> shape_second_derivative(Precision r, Precision s) const override {
        StaticMatrix<3, 3> ddN;
        ddN << 0, 0, 0,
               0, 0, 0,
               0, 0, 0;
        return ddN;
    }

    StaticMatrix<3, 2> node_coords_local() const override {
        StaticMatrix<3, 2> local_coords;
        local_coords << 0.0, 0.0,   // Node 1
                        1.0, 0.0,   // Node 2
                        0.0, 1.0;   // Node 3
        return local_coords;
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

    StaticVector<2> closest_point_on_boundary(const StaticVector<3>& global, const StaticMatrix<3, 3>& node_coords) const {
        Line2B line1({0, 1});
        Line2B line2({1, 2});
        Line2B line3({2, 0});

        Precision line1_p = line1.global_to_local(global, node_coords);
        Precision line2_p = line2.global_to_local(global, node_coords);
        Precision line3_p = line3.global_to_local(global, node_coords);

        StaticVector<3> p1 = line1.local_to_global(line1_p, node_coords);
        StaticVector<3> p2 = line2.local_to_global(line2_p, node_coords);
        StaticVector<3> p3 = line3.local_to_global(line3_p, node_coords);

        Precision d1 = (p1 - global).squaredNorm();
        Precision d2 = (p2 - global).squaredNorm();
        Precision d3 = (p3 - global).squaredNorm();

        if (d1 < d2 && d1 < d3) return {line1_p, 0};
        if (d3 < d1 && d3 < d2) return {0, line3_p};
        return {1-line2_p, line2_p};
    }

    bool in_bounds(const StaticVector<2>& local) const {
        return local(0) >= 0 && local(1) >= 0 && local(0) + local(1) <= 1;
    }

};

}  // namespace fem::model
