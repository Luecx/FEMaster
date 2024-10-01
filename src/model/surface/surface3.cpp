/******************************************************************************
 * @file surface3.cpp
 * @brief Implementation file for the Surface3 class.
 *
 * @details Contains the definitions of member functions for the `Surface3` class,
 *          which represents a triangular surface element with three nodes. This
 *          file includes methods for shape function evaluation, shape function
 *          derivatives, and element boundary projections.
 *
 * @date Created on 27.09.2024 by Finn Eggers
 ******************************************************************************/

#include "surface3.h"

/**
 * @brief Constructor for the Surface3 class (triangular surface element with 3 nodes).
 *
 * @param pNodeIds Array of node IDs corresponding to the surface element nodes.
 */
fem::model::Surface3::Surface3(const std::array<ID, 3>& pNodeIds)
    : SurfaceInterface<3>(pNodeIds) {}

/**
 * @brief Compute the shape functions for the triangular element.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<3, 1> Vector of shape function values.
 */
StaticMatrix<3, 1> fem::model::Surface3::shape_function(Precision r, Precision s) const {
    StaticMatrix<3, 1> N;
    N(0, 0) = 1 - r - s;    // Shape function for node 1
    N(1, 0) = r;            // Shape function for node 2
    N(2, 0) = s;            // Shape function for node 3
    return N;
}

/**
 * @brief Compute the first derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<3, 2> Matrix of shape function derivatives.
 */
StaticMatrix<3, 2> fem::model::Surface3::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<3, 2> dN;
    dN(0, 0) = -1;
    dN(0, 1) = -1;    // dN1/dr, dN1/ds
    dN(1, 0) = 1;
    dN(1, 1) = 0;    // dN2/dr, dN2/ds
    dN(2, 0) = 0;
    dN(2, 1) = 1;    // dN3/dr, dN3/ds
    return dN;
}

/**
 * @brief Compute the second derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<3, 3> Zero matrix for second-order derivatives.
 */
StaticMatrix<3, 3> fem::model::Surface3::shape_second_derivative(Precision r, Precision s) const {
    StaticMatrix<3, 3> ddN;
    ddN << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    return ddN;
}

/**
 * @brief Retrieve the local coordinates of the nodes for the triangular element.
 *
 * @return StaticMatrix<3, 2> Matrix of local node coordinates.
 */
StaticMatrix<3, 2> fem::model::Surface3::node_coords_local() const {
    StaticMatrix<3, 2> local_coords;
    local_coords << 0.0, 0.0,    // Node 1
                    1.0, 0.0,    // Node 2
                    0.0, 1.0;    // Node 3
    return local_coords;
}

/**
 * @brief Return the integration scheme for the element.
 *
 * @return const quadrature::Quadrature& Quadrature scheme for the triangular element.
 */
const fem::quadrature::Quadrature& fem::model::Surface3::integration_scheme() const {
    static const quadrature::Quadrature quad {quadrature::DOMAIN_ISO_TRI, quadrature::ORDER_LINEAR};
    return quad;
}

/**
 * @brief Compute the closest point on the boundary to a given global point.
 *
 * @param global Global coordinates of the point.
 * @param node_coords Coordinates of the element nodes.
 * @return Vec2 Local coordinates of the closest boundary point.
 */
Vec2 fem::model::Surface3::closest_point_on_boundary(const StaticVector<3>& global, const StaticMatrix<3, 3>& node_coords) const {
    // Implement boundary projections using line elements
    // Line elements are defined between nodes
    Line2B line1({0, 1});
    Line2B line2({1, 2});
    Line2B line3({2, 0});

    // Compute projections onto the three boundary lines
    Precision line1_p = line1.global_to_local(global, node_coords);
    Precision line2_p = line2.global_to_local(global, node_coords);
    Precision line3_p = line3.global_to_local(global, node_coords);

    // Convert local line parameters to global points
    StaticVector<3> p1 = line1.local_to_global(line1_p, node_coords);
    StaticVector<3> p2 = line2.local_to_global(line2_p, node_coords);
    StaticVector<3> p3 = line3.local_to_global(line3_p, node_coords);

    // Calculate squared distances
    Precision d1 = (p1 - global).squaredNorm();
    Precision d2 = (p2 - global).squaredNorm();
    Precision d3 = (p3 - global).squaredNorm();

    // Return the local coordinates of the closest boundary point
    if (d1 <= d2 && d1 <= d3) return {line1_p, 0};
    if (d3 <= d1 && d3 <= d2) return {0, 1-line3_p};
    return {1 - line2_p, line2_p};
}

/**
 * @brief Check if a local point is within the bounds of the triangular element.
 *
 * @param local Local coordinates (r, s).
 * @return bool True if the point is within bounds, false otherwise.
 */
bool fem::model::Surface3::in_bounds(const Vec2& local) const {
    return local(0) >= 0 && local(1) >= 0 && local(0) + local(1) <= 1;
}
