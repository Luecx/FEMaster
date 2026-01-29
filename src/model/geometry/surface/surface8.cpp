/**
 * @file surface8.cpp
 * @brief Implements the functions for the Surface8 class.
 *
 * @details This file contains the implementation of shape functions, shape
 *          derivatives, local to global coordinate mappings, and surface
 *          integration schemes for the 8-node quadrilateral surface element.
 *
 * @date Created on 01.10.2024 by Finn Eggers
 */

#include "surface8.h"
#include "../line/line3a.h"

namespace fem::model {

/**
 * @brief Constructor for the Surface8 class (quadrilateral surface element with 8 nodes).
 *
 * @param pNodeIds Array of node IDs corresponding to the surface element nodes.
 */
fem::model::Surface8::Surface8(const std::array<ID, 8>& pNodeIds)
    : Surface<8>(pNodeIds) {}

/**
 * @brief Compute the shape functions for the 8-node quadrilateral element.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<8, 1> Vector of shape function values.
 */
StaticMatrix<8, 1> Surface8::shape_function(Precision r, Precision s) const {
    StaticMatrix<8, 1> N;
    N(0, 0) = 0.25 * (1 - r) * (1 - s) * (-1 - r - s);   // Node 1
    N(1, 0) = 0.25 * (1 + r) * (1 - s) * (-1 + r - s);   // Node 2
    N(2, 0) = 0.25 * (1 + r) * (1 + s) * (-1 + r + s);   // Node 3
    N(3, 0) = 0.25 * (1 - r) * (1 + s) * (-1 - r + s);   // Node 4
    N(4, 0) = 0.5 * (1 - r*r) * (1 - s);                 // Node 5 (mid-edge)
    N(5, 0) = 0.5 * (1 + r) * (1 - s*s);                 // Node 6 (mid-edge)
    N(6, 0) = 0.5 * (1 - r*r) * (1 + s);                 // Node 7 (mid-edge)
    N(7, 0) = 0.5 * (1 - r) * (1 - s*s);                 // Node 8 (mid-edge)
    return N;
}

/**
 * @brief Compute the first derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<8, 2> Matrix of shape function derivatives.
 */
StaticMatrix<8, 2> Surface8::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<8, 2> dN;

    // Derivatives of the shape functions with respect to r and s
    dN(0, 0) =  0.25 * (-2*r - s) * (s - 1);  dN(0, 1) =  0.25 * (-r - 2*s) * (r - 1);  // dN1/dr, dN1/ds
    dN(1, 0) =  0.25 * (-2*r + s) * (s - 1);  dN(1, 1) =  0.25 * (-r + 2*s) * (r + 1);  // dN2/dr, dN2/ds
    dN(2, 0) =  0.25 * ( 2*r + s) * (s + 1);  dN(2, 1) =  0.25 * (r + 1) * (r + 2*s);   // dN3/dr, dN3/ds
    dN(3, 0) =  0.25 * ( 2*r - s) * (s + 1);  dN(3, 1) =  0.25 * (r - 1) * (r - 2*s);   // dN4/dr, dN4/ds
    dN(4, 0) =  r * (s - 1);                  dN(4, 1) =  0.5 * (r*r - 1);              // dN5/dr, dN5/ds
    dN(5, 0) =  0.5 * (1 - s*s);              dN(5, 1) = -1.0 * s * (r + 1);            // dN6/dr, dN6/ds
    dN(6, 0) = -r * (1 + s);                  dN(6, 1) =  0.5 * (1 - r*r);              // dN7/dr, dN7/ds
    dN(7, 0) =  0.5 * (s*s - 1);              dN(7, 1) =  1.0 * s * (r - 1);            // dN8/dr, dN8/ds


    return dN;
}

/**
 * @brief Compute the second derivatives of the shape functions.
 *
 * @param r Local coordinate in the parametric space.
 * @param s Local coordinate in the parametric space.
 * @return StaticMatrix<8, 3> Matrix of second-order shape function derivatives.
 */
StaticMatrix<8, 3> Surface8::shape_second_derivative(Precision r, Precision s) const {
    StaticMatrix<8, 3> ddN;

    // Second derivatives of the shape functions
    ddN(0, 0) =  0.5 - 0.5 * s; ddN(0, 1) =  0.5 - 0.5 * r; ddN(0, 2) = -0.5 * r - 0.5 * s + 0.25;  // d²N1/dr², d²N1/ds², d²N1/(drds)
    ddN(1, 0) =  0.5 - 0.5 * s; ddN(1, 1) =  0.5 * r + 0.5; ddN(1, 2) = -0.5 * r + 0.5 * s - 0.25;  // d²N2/dr², d²N2/ds², d²N2/(drds)
    ddN(2, 0) =  0.5 * s + 0.5; ddN(2, 1) =  0.5 * r + 0.5; ddN(2, 2) =  0.5 * r + 0.5 * s + 0.25;  // d²N3/dr², d²N3/ds², d²N3/(drds)
    ddN(3, 0) =  0.5 * s + 0.5; ddN(3, 1) =  0.5 - 0.5 * r; ddN(3, 2) =  0.5 * r - 0.5 * s - 0.25;  // d²N4/dr², d²N4/ds², d²N4/(drds)
    ddN(4, 0) =  1.0 * s - 1.0; ddN(4, 1) =  0.0;           ddN(4, 2) =  1.0 * r;                   // d²N5/dr², d²N5/ds², d²N5/(drds)
    ddN(5, 0) =  0.0;           ddN(5, 1) = -1.0 * r - 1.0; ddN(5, 2) = -1.0 * s;                   // d²N6/dr², d²N6/ds², d²N6/(drds)
    ddN(6, 0) = -1.0 * s - 1.0; ddN(6, 1) =  0.0;           ddN(6, 2) = -1.0 * r;                   // d²N7/dr², d²N7/ds², d²N7/(drds)
    ddN(7, 0) =  0.0;           ddN(7, 1) =  1.0 * r - 1.0; ddN(7, 2) =  1.0 * s;                   // d²N8/dr², d²N8/ds², d²N8/(drds)

    return ddN;
}

/**
 * @brief Retrieve the local coordinates of the nodes for the 8-node quadrilateral element.
 *
 * @return StaticMatrix<8, 2> Local coordinates of the nodes.
 */
StaticMatrix<8, 2> Surface8::node_coords_local() const {
    StaticMatrix<8, 2> local_coords;
    local_coords << -1.0, -1.0,   // Node 1
                     1.0, -1.0,   // Node 2
                     1.0,  1.0,   // Node 3
                    -1.0,  1.0,   // Node 4
                     0.0, -1.0,   // Node 5 (mid-edge)
                     1.0,  0.0,   // Node 6 (mid-edge)
                     0.0,  1.0,   // Node 7 (mid-edge)
                    -1.0,  0.0;   // Node 8 (mid-edge)
    return local_coords;
}

/**
 * @brief Return the integration scheme for the 8-node quadrilateral element.
 *
 * @return const quadrature::Quadrature& Quadrature scheme for the 8-node element.
 */
const fem::quadrature::Quadrature& Surface8::integration_scheme() const {
    static const quadrature::Quadrature quad{quadrature::DOMAIN_ISO_QUAD, quadrature::ORDER_QUADRATIC};
    return quad;
}

/**
 * @brief Compute the closest point on the boundary to a given global point.
 *
 * @param global Global coordinates of the point.
 * @param node_coords Coordinates of the element nodes.
 * @return Vec2 Local coordinates of the closest boundary point.
 */
Vec2 Surface8::closest_point_on_boundary(const Vec3& global, const StaticMatrix<8, 3>& node_coords) const {
    // Boundary checks using line elements defined between nodes
    Line3A line1({0, 1, 4});  // Line from node 1 to node 2
    Line3A line2({1, 2, 5});  // Line from node 2 to node 3
    Line3A line3({2, 3, 6});  // Line from node 3 to node 4
    Line3A line4({3, 0, 7});  // Line from node 4 to node 1

    Field node_field("SURFACE8_BOUNDARY", FieldDomain::NODE, 8, 3);
    for (Index i = 0; i < 8; ++i) {
        for (Index j = 0; j < 3; ++j) {
            node_field(i, j) = node_coords(i, j);
        }
    }

    // Compute projections onto the four boundary lines
    Precision line1_p = line1.global_to_local(global, node_field);
    Precision line2_p = line2.global_to_local(global, node_field);
    Precision line3_p = line3.global_to_local(global, node_field);
    Precision line4_p = line4.global_to_local(global, node_field);

    // Convert local line parameters to global points
    Vec3 p1 = line1.local_to_global(line1_p, node_field);
    Vec3 p2 = line2.local_to_global(line2_p, node_field);
    Vec3 p3 = line3.local_to_global(line3_p, node_field);
    Vec3 p4 = line4.local_to_global(line4_p, node_field);

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
bool Surface8::in_bounds(const Vec2& local) const {
    return local(0) >= -1 && local(0) <= 1 && local(1) >= -1 && local(1) <= 1;
}

}  // namespace fem::model
