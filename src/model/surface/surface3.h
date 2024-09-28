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

    // Function signature for the method in your class
    StaticVector<2> global_to_local(const StaticVector<3>& global, const NodeData& node_coords_system, bool clip) const override {
        // Local coordinate vector to return
        StaticVector<2> local;

        // Retrieve node coordinates (adjust node_coords row selection based on your class)
        auto node_coords = this->node_coords_global(node_coords_system);

        // Define vectors spanning the plane from the first node
        StaticVector<3> v0 = node_coords.row(1) - node_coords.row(0);  // Node 2 - Node 1
        StaticVector<3> v1 = node_coords.row(2) - node_coords.row(0);  // Node 3 - Node 1
        StaticVector<3> v2 = global - node_coords.row(0).transpose();  // From Node 1 to the point

        // Step 1: Project `v2` onto the plane defined by `v0` and `v1`
        // Calculate the normal vector to the plane
        StaticVector<3> normal = v0.cross(v1).normalized();

        // Project `v2` onto the plane: projected_vector = v2 - (v2 • normal) * normal
        StaticVector<3> v2_proj = v2 - (v2.dot(normal)) * normal;

        // Step 2: Convert to local coordinates in the plane spanned by `v0` and `v1` using Normal Equations
        // Set up matrix U = [v0 | v1]
        StaticMatrix<3, 2> U;
        U.col(0) = v0;
        U.col(1) = v1;

        // Calculate U^T * U
        Eigen::Matrix<Precision, 2, 2> UtU = U.transpose() * U;

        // Calculate U^T * v2_proj
        StaticVector<2> Utv2_proj = U.transpose() * v2_proj;

        // Step 3: Use matrix inversion to solve the equation: local = (U^T * U)^(-1) * (U^T * v2_proj)
        local = UtU.inverse() * Utv2_proj;

        // Step 4: Clip the coordinates if needed (if the point is outside the triangle)
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
