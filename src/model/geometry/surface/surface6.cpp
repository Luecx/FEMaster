/**
 * @file surface6.cpp
 * @brief Implements the functions for the Surface6 class.
 *
 * @details This file contains the implementation of shape functions, shape
 *          derivatives, local to global coordinate mappings, and surface
 *          integration schemes for the 6-node triangular surface element.
 *
 * @date Created on 01.10.2024 by Finn Eggers
 */

#include "surface6.h"
#include "../line/line3b.h"

namespace fem::model {

/**
 * @brief Constructor for the Surface6 class (quadrilateral surface element with 6 nodes).
 *
 * @param pNodeIds Array of node IDs corresponding to the surface element nodes.
 */
fem::model::Surface6::Surface6(const std::array<ID, 6>& pNodeIds)
    : Surface<6>(pNodeIds) {}

/**
 * @brief Compute the shape functions for the 6-node triangular element.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<6, 1> Vector of shape function values.
 */
StaticMatrix<6, 1> Surface6::shape_function(Precision r, Precision s) const {
    StaticMatrix<6, 1> N;
    N(0, 0) = 1 - 3*(r + s) + 2*(r + s)*(r + s);  // Shape function for node 1
    N(1, 0) = r * (2*r - 1);                       // Shape function for node 2
    N(2, 0) = s * (2*s - 1);                       // Shape function for node 3
    N(3, 0) = 4 * r * (1 - r - s);                 // Shape function for node 4
    N(4, 0) = 4 * r * s;                           // Shape function for node 5
    N(5, 0) = 4 * s * (1 - r - s);                 // Shape function for node 6
    return N;
}

/**
 * @brief Compute the first derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<6, 2> Matrix of shape function derivatives.
 */
StaticMatrix<6, 2> Surface6::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<6, 2> dN;
    dN(0, 0) = -3 + 4*(r + s); dN(0, 1) = -3 + 4*(r + s);
    dN(1, 0) = 4*r - 1;        dN(1, 1) = 0;
    dN(2, 0) = 0;              dN(2, 1) = 4*s - 1;
    dN(3, 0) = 4 - 8*r - 4*s;  dN(3, 1) = -4*r;
    dN(4, 0) = 4*s;            dN(4, 1) = 4*r;
    dN(5, 0) = -4*s;           dN(5, 1) = 4 - 4*r - 8*s;
    return dN;
}

/**
 * @brief Compute the second derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<6, 3> Matrix of second-order shape function derivatives.
 */
StaticMatrix<6, 3> Surface6::shape_second_derivative(Precision r, Precision s) const {
    (void) r;
    (void) s;

    StaticMatrix<6, 3> ddN;

    // N1 = 1 - 3(r + s) + 2(r + s)^2
    ddN(0, 0) = 4;              // ∂²N1/∂r²
    ddN(0, 1) = 4;              // ∂²N1/∂s²
    ddN(0, 2) = 4;              // ∂²N1/∂(rs)

    // N2 = r(2r - 1)
    ddN(1, 0) = 4;              // ∂²N2/∂r²
    ddN(1, 1) = 0;              // ∂²N2/∂s²
    ddN(1, 2) = 0;              // ∂²N2/∂(rs)

    // N3 = s(2s - 1)
    ddN(2, 0) = 0;              // ∂²N3/∂r²
    ddN(2, 1) = 4;              // ∂²N3/∂s²
    ddN(2, 2) = 0;              // ∂²N3/∂(rs)

    // N4 = 4r(1 - r - s)
    ddN(3, 0) = -8;             // ∂²N4/∂r²
    ddN(3, 1) = 0;              // ∂²N4/∂s²
    ddN(3, 2) = -4;             // ∂²N4/∂(rs)

    // N5 = 4rs
    ddN(4, 0) = 0;              // ∂²N5/∂r²
    ddN(4, 1) = 0;              // ∂²N5/∂s²
    ddN(4, 2) = 4;              // ∂²N5/∂(rs)

    // N6 = 4s(1 - r - s)
    ddN(5, 0) = 0;              // ∂²N6/∂r²
    ddN(5, 1) = -8;             // ∂²N6/∂s²
    ddN(5, 2) = -4;             // ∂²N6/∂(rs)

    return ddN;
}
/**
 * @brief Retrieve the local coordinates of the nodes for the 6-node triangular element.
 *
 * @return StaticMatrix<6, 2> Local coordinates of the nodes.
 */
StaticMatrix<6, 2> Surface6::node_coords_local() const {
    StaticMatrix<6, 2> local_coords;
    local_coords << 0, 0,
                    1, 0,
                    0, 1,
                    0.5, 0,
                    0.5, 0.5,
                    0, 0.5;
    return local_coords;
}

/**
 * @brief Return the integration scheme for the element.
 *
 * @return const quadrature::Quadrature& Quadrature scheme for the 6-node triangular element.
 */
const fem::quadrature::Quadrature& Surface6::integration_scheme() const {
    static const quadrature::Quadrature quad{quadrature::DOMAIN_ISO_TRI, quadrature::ORDER_QUADRATIC};
    return quad;
}

/**
 * @brief Compute the closest point on the boundary to a given global point.
 *
 * @param global Global coordinates of the point.
 * @param node_coords Coordinates of the element nodes.
 * @return Vec2 Local coordinates of the closest boundary point.
 */
Vec2 Surface6::closest_point_on_boundary(const Vec3& global, const StaticMatrix<6, 3>& node_coords) const {
    // Implement boundary projections using line elements
    // Line elements are defined between nodes
    Line3B line1({0, 1, 3});
    Line3B line2({1, 2, 4});
    Line3B line3({2, 0, 5});

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
bool Surface6::in_bounds(const Vec2& local) const {
    return local(0) >= 0 && local(1) >= 0 && local(0) + local(1) <= 1;
}

}  // namespace fem::model
