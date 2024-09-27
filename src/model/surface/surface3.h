#pragma once

#include "surface.h"

#include <array>

namespace fem::model {

/******************************************************************************
 * @class Surface3
 * @brief Derived class for a triangular surface element with 3 nodes (S3D3).
 ******************************************************************************/

struct Surface3 : public SurfaceInterface<3> {

    Surface3(const std::array<ID, 3>& pNodeIds)
        : SurfaceInterface<3>(pNodeIds) {}

    /******************************************************************************
     * @brief Computes the shape function for a triangular element with 3 nodes.
     *        Uses natural coordinates r and s.
     *
     * @param r The first natural coordinate (ranges from 0 to 1).
     * @param s The second natural coordinate (ranges from 0 to 1).
     * @return Shape function values evaluated at (r, s).
     ******************************************************************************/
    StaticMatrix<3, 1> shape_function(Precision r, Precision s) const override {
        StaticMatrix<3, 1> N;
        N(0, 0) = 1 - r - s;    // Shape function for node 1
        N(1, 0) = r;            // Shape function for node 2
        N(2, 0) = s;            // Shape function for node 3
        return N;
    }

    /******************************************************************************
     * @brief Computes the shape function derivatives with respect to r and s.
     *
     * @param r The first natural coordinate.
     * @param s The second natural coordinate.
     * @return A matrix of shape derivatives.
     ******************************************************************************/
    StaticMatrix<3, 2> shape_derivative(Precision r, Precision s) const override {
        StaticMatrix<3, 2> dN;
        dN(0, 0) = -1;
        dN(0, 1) = -1;    // dN1/dr, dN1/ds
        dN(1, 0) = 1;
        dN(1, 1) = 0;    // dN2/dr, dN2/ds
        dN(2, 0) = 0;
        dN(2, 1) = 1;    // dN3/dr, dN3/ds
        return dN;
    }

    /******************************************************************************
     * @brief Gets the local coordinates of the triangle's nodes.
     *
     * @return A matrix containing the local coordinates of the 3 nodes.
     ******************************************************************************/
    StaticMatrix<3, 2> node_coords_local() const override {
        StaticMatrix<3, 2> local_coords;
        local_coords << 0.0, 0.0,    // Node 1 (corner)
            1.0, 0.0,                // Node 2 (corner)
            0.0, 1.0;                // Node 3 (corner)
        return local_coords;
    }

    /******************************************************************************
     * @brief Maps a global coordinate to the local coordinate system.
     *
     * @param global The global coordinate to be mapped.
     * @return The corresponding local coordinates (r, s).
     ******************************************************************************/
    StaticVector<2>
        global_to_local(const StaticVector<3>& global, const NodeData& node_coords_system, bool clip) const override {
        // For a linear triangle, we can use an analytical solution for projection.
        StaticVector<2> local;
        auto            node_coords = this->node_coords_global(node_coords_system);

        // Compute vectors for the triangle's edges
        StaticVector<3> v0 = node_coords.row(2) - node_coords.row(0);    // Edge 1 (Node 3 - Node 1)
        StaticVector<3> v1 = node_coords.row(1) - node_coords.row(0);    // Edge 2 (Node 2 - Node 1)
        StaticVector<3> v2 = global - node_coords.row(0).transpose();    // Vector from Node 1 to the given point

        // Compute coefficients to solve the local coordinates
        Precision d00   = v0.dot(v0);
        Precision d01   = v0.dot(v1);
        Precision d11   = v1.dot(v1);
        Precision d20   = v2.dot(v0);
        Precision d21   = v2.dot(v1);
        Precision denom = d00 * d11 - d01 * d01;

        // Calculate barycentric coordinates (r, s)
        local(0) = (d00 * d21 - d01 * d20) / denom;    // r-coordinate
        local(1) = (d11 * d20 - d01 * d21) / denom;    // s-coordinate

        // Clip the coordinates if requested
        if (clip) {
            // Check if the point is inside the triangle
            if (local(0) < 0 || local(1) < 0 || local(0) + local(1) > 1) {
                // Project the point onto the edges of the triangle and find the closest point
                StaticVector<2> clipped_local = local;

                // Node positions in local coordinates
                StaticVector<2> p0(0, 0);    // Node 1
                StaticVector<2> p1(1, 0);    // Node 2
                StaticVector<2> p2(0, 1);    // Node 3

                // Project onto the three edges and find the closest point
                clipped_local = closest_point_on_triangle_edges(local, p0, p1, p2);

                return clipped_local;
            }
        }

        return local;
    }

    /******************************************************************************
     * @brief Computes the closest point on a triangle's edges in barycentric coordinates.
     *
     * @param local Barycentric coordinates to be clipped.
     * @param p0 Local coordinates of the first node.
     * @param p1 Local coordinates of the second node.
     * @param p2 Local coordinates of the third node.
     * @return Clipped barycentric coordinates.
     ******************************************************************************/
    StaticVector<2> closest_point_on_triangle_edges(const StaticVector<2>& local,
                                                    const StaticVector<2>& p0,
                                                    const StaticVector<2>& p1,
                                                    const StaticVector<2>& p2) const {
        // Calculate the closest point on each of the triangle edges
        StaticVector<2> c0;
        StaticVector<2> c1;
        StaticVector<2> c2;

        bool            outside0, outside1, outside2;

        std::tie(c0, outside0) = closest_point_on_edge(local, p0, p1);    // Closest on Edge 1 (p0-p1)
        std::tie(c1, outside1) = closest_point_on_edge(local, p1, p2);    // Closest on Edge 2 (p1-p2)
        std::tie(c2, outside2) = closest_point_on_edge(local, p2, p0);    // Closest on Edge 3 (p2-p0)

        // Compute the squared distances to the original point
        Precision d0 = (local - c0).squaredNorm();
        Precision d1 = (local - c1).squaredNorm();
        Precision d2 = (local - c2).squaredNorm();

        // If the point is outside all edges, use the vertices as the closest points
        if (outside0 && outside1 && outside2) {
            Precision d_p0 = (local - p0).squaredNorm();
            Precision d_p1 = (local - p1).squaredNorm();
            Precision d_p2 = (local - p2).squaredNorm();

            if (d_p0 < d_p1 && d_p0 < d_p2) {
                return p0;    // Closest to vertex 1
            } else if (d_p1 < d_p2) {
                return p1;    // Closest to vertex 2
            } else {
                return p2;    // Closest to vertex 3
            }
        }

        // Return the closest point on the triangle edges
        if (d0 < d1 && d0 < d2) {
            return c0;
        } else if (d1 < d2) {
            return c1;
        } else {
            return c2;
        }
    }

    /******************************************************************************
     * @brief Finds the closest point on a line segment defined by two endpoints.
     *
     * @param local Barycentric coordinates to be clipped.
     * @param p0 Local coordinates of the first endpoint.
     * @param p1 Local coordinates of the second endpoint.
     * @return Pair of clipped barycentric coordinates along the edge and a boolean flag.
     *         If the flag is true, the closest point is outside the segment.
     ******************************************************************************/
    std::pair<StaticVector<2>, bool> closest_point_on_edge(const StaticVector<2>& local,
                                                           const StaticVector<2>& p0,
                                                           const StaticVector<2>& p1) const {
        // Vector from p0 to p1
        StaticVector<2> edge = p1 - p0;

        // Compute the projection of the point onto the edge
        Precision t = (local - p0).dot(edge) / edge.squaredNorm();

        // Determine if the projection is outside the edge segment
        bool outside = t < 0.0 || t > 1.0;

        // Clamp t to [0, 1] to stay on the edge segment
        t = std::max(0.0, std::min(1.0, t));

        // Compute the closest point along the edge
        StaticVector<2> closest_point = p0 + t * edge;

        return std::make_pair(closest_point, outside);
    }

    /******************************************************************************
     * @brief Returns the integration scheme for a 3-node triangle.
     *
     * @return The quadrature scheme to be used for integration.
     ******************************************************************************/
    const quadrature::Quadrature& integration_scheme() const override {
        const static quadrature::Quadrature quad {quadrature::DOMAIN_ISO_TRI, quadrature::ORDER_LINEAR};
        return quad;
    }
};

}    // namespace fem::model
