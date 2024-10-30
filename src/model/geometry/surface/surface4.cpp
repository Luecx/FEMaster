/******************************************************************************
 * @file surface4.cpp
 * @brief Implementation file for the Surface4 class.
 *
 * @details Contains the definitions of member functions for the `Surface4` class,
 *          which represents a quadrilateral surface element with four nodes. This
 *          file includes methods for shape function evaluation, shape function
 *          derivatives, and element boundary projections.
 *
 * @date Created on 27.09.2024 by Finn Eggers
 ******************************************************************************/

#include "surface4.h"

namespace fem::model {
/**
 * @brief Constructor for the Surface4 class (quadrilateral surface element with 4 nodes).
 *
 * @param pNodeIds Array of node IDs corresponding to the surface element nodes.
 */
Surface4::Surface4(const std::array<ID, 4>& pNodeIds)
    : Surface<4>(pNodeIds) {}

/**
 * @brief Compute the shape functions for the quadrilateral element.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<4, 1> Vector of shape function values.
 */
StaticMatrix<4, 1> Surface4::shape_function(Precision r, Precision s) const {
    StaticMatrix<4, 1> N;
    N(0, 0) = 0.25 * (1 - r) * (1 - s);  // Shape function for node 1
    N(1, 0) = 0.25 * (1 + r) * (1 - s);  // Shape function for node 2
    N(2, 0) = 0.25 * (1 + r) * (1 + s);  // Shape function for node 3
    N(3, 0) = 0.25 * (1 - r) * (1 + s);  // Shape function for node 4
    return N;
}

/**
 * @brief Compute the first derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<4, 2> Matrix of shape function derivatives.
 */
StaticMatrix<4, 2> Surface4::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<4, 2> dN;
    dN(0, 0) = -0.25 * (1 - s);  dN(0, 1) = -0.25 * (1 - r);  // dN1/dr, dN1/ds
    dN(1, 0) =  0.25 * (1 - s);  dN(1, 1) = -0.25 * (1 + r);  // dN2/dr, dN2/ds
    dN(2, 0) =  0.25 * (1 + s);  dN(2, 1) =  0.25 * (1 + r);  // dN3/dr, dN3/ds
    dN(3, 0) = -0.25 * (1 + s);  dN(3, 1) =  0.25 * (1 - r);  // dN4/dr, dN4/ds
    return dN;
}

/**
 * @brief Compute the second derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<4, 3> Matrix of second-order shape function derivatives.
 */
StaticMatrix<4, 3> Surface4::shape_second_derivative(Precision r, Precision s) const {
    (void) r;
    (void) s;
    StaticMatrix<4, 3> ddN;
    ddN << 0, 0,  0.25,  // d²N1/dr², d²N1/ds², d²N1/(drds)
           0, 0, -0.25,  // d²N2/dr², d²N2/ds², d²N2/(drds)
           0, 0,  0.25,  // d²N3/dr², d²N3/ds², d²N3/(drds)
           0, 0, -0.25;  // d²N4/dr², d²N4/ds², d²N4/(drds)
    return ddN;
}

/**
 * @brief Retrieve the local coordinates of the nodes for the quadrilateral element.
 *
 * @return StaticMatrix<4, 2> Local coordinates of the nodes.
 */
StaticMatrix<4, 2> Surface4::node_coords_local() const {
    StaticMatrix<4, 2> local_coords;
    local_coords << -1.0, -1.0,  // Node 1
                     1.0, -1.0,  // Node 2
                     1.0,  1.0,  // Node 3
                    -1.0,  1.0;  // Node 4
    return local_coords;
}

/**
 * @brief Return the integration scheme for the element.
 *
 * @return const quadrature::Quadrature& Quadrature scheme for the quadrilateral element.
 */
const fem::quadrature::Quadrature& Surface4::integration_scheme() const {
    static const quadrature::Quadrature quad{quadrature::DOMAIN_ISO_QUAD, quadrature::ORDER_LINEAR};
    return quad;
}

/**
 * @brief Compute the closest point on the boundary to a given global point.
 *
 * @param global Global coordinates of the point.
 * @param node_coords Coordinates of the element nodes.
 * @return Vec2 Local coordinates of the closest boundary point.
 */
Vec2 Surface4::closest_point_on_boundary(const Vec3& global, const StaticMatrix<4, 3>& node_coords) const {
    // Boundary checks using line elements defined between nodes
    Line2A line1({0, 1});  // Line from node 1 to node 2
    Line2A line2({1, 2});  // Line from node 2 to node 3
    Line2A line3({2, 3});  // Line from node 3 to node 4
    Line2A line4({3, 0});  // Line from node 4 to node 1

    // Compute projections onto the four boundary lines
    Precision line1_p = line1.global_to_local(global, node_coords);
    Precision line2_p = line2.global_to_local(global, node_coords);
    Precision line3_p = line3.global_to_local(global, node_coords);
    Precision line4_p = line4.global_to_local(global, node_coords);

    // Convert local line parameters to global points
    Vec3 p1 = line1.local_to_global(line1_p, node_coords);
    Vec3 p2 = line2.local_to_global(line2_p, node_coords);
    Vec3 p3 = line3.local_to_global(line3_p, node_coords);
    Vec3 p4 = line4.local_to_global(line4_p, node_coords);

    // Calculate squared distances
    Precision d1 = (p1 - global).squaredNorm();
    Precision d2 = (p2 - global).squaredNorm();
    Precision d3 = (p3 - global).squaredNorm();
    Precision d4 = (p4 - global).squaredNorm();

    // Return the local coordinates of the closest boundary point
    if (d1 <= d2 && d1 <= d3 && d1 <= d4) return {line1_p, -1};
    if (d2 <= d1 && d2 <= d3 && d2 <= d4) return {1, line2_p};
    if (d3 <= d1 && d3 <= d2 && d3 <= d4) return {-line3_p, 1};
    return {-1, -line4_p};
}

/**
 * @brief Check if a local point is within the bounds of the quadrilateral element.
 *
 * @param local Local coordinates (r, s).
 * @return bool True if the point is within bounds, false otherwise.
 */
bool Surface4::in_bounds(const Vec2& local) const {
    return local(0) >= -1 && local(0) <= 1 && local(1) >= -1 && local(1) <= 1;
}

}  // namespace fem::model

