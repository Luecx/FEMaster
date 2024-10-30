//
// Created by Finn Eggers on 27.08.2024.
//

#include "c3d4.h"

namespace fem {
namespace model {

//-----------------------------------------------------------------------------
// C3D4 Constructor
//-----------------------------------------------------------------------------
/**
 * @brief Constructor that initializes the C3D4 element with an element ID
 * and the corresponding node IDs.
 *
 * @param pElemId The unique ID of the element.
 * @param pNodeIds Array containing IDs of the 4 nodes that define the element.
 */
C3D4::C3D4(ID pElemId, const std::array<ID, 4>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

//-----------------------------------------------------------------------------
// Shape Function
//-----------------------------------------------------------------------------
/**
 * @brief Computes the shape functions for the 4-node tetrahedral element at
 * the given local coordinates (r, s, t).
 *
 * The shape functions are defined as follows:
 * - N1 = r
 * - N2 = s
 * - N3 = 1 - r - s - t
 * - N4 = t
 *
 * @param r Local coordinate along the r-direction.
 * @param s Local coordinate along the s-direction.
 * @param t Local coordinate along the t-direction.
 * @return StaticMatrix<4, 1> The computed shape function values at the
 * specified local coordinates.
 */
StaticMatrix<4, 1> C3D4::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<4, 1>  res {};

    // Define the shape functions for the 4 nodes
    res(0) = r;  // Shape function N1
    res(1) = s;  // Shape function N2
    res(2) = 1 - r - s - t;  // Shape function N3
    res(3) = t;  // Shape function N4

    return res;
}

//-----------------------------------------------------------------------------
// Shape Function Derivatives
//-----------------------------------------------------------------------------
/**
 * @brief Computes the derivatives of the shape functions with respect to the
 * local coordinates (r, s, t).
 *
 * The derivatives of the shape functions are given by:
 * - dN1/dr = 1, dN1/ds = 0, dN1/dt = 0
 * - dN2/dr = 0, dN2/ds = 1, dN2/dt = 0
 * - dN3/dr = -1, dN3/ds = -1, dN3/dt = -1
 * - dN4/dr = 0, dN4/ds = 0, dN4/dt = 1
 *
 * @param r Local coordinate along the r-direction.
 * @param s Local coordinate along the s-direction.
 * @param t Local coordinate along the t-direction.
 * @return StaticMatrix<4, 3> The derivatives of the shape functions
 * evaluated at the local coordinates.
 */
StaticMatrix<4, 3> C3D4::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<4, 3> local_shape_derivative {};
    local_shape_derivative.setZero();

    // Define the derivatives of the shape functions
    local_shape_derivative(0, 0) = 1;  // dN1/dr

    local_shape_derivative(1, 1) = 1;  // dN2/ds

    local_shape_derivative(2, 0) = -1; // dN3/dr
    local_shape_derivative(2, 1) = -1; // dN3/ds
    local_shape_derivative(2, 2) = -1; // dN3/dt

    local_shape_derivative(3, 2) = 1;  // dN4/dt

    return local_shape_derivative;
}

//-----------------------------------------------------------------------------
// Node Local Coordinates
//-----------------------------------------------------------------------------
/**
 * @brief Provides the local coordinates of the nodes for the C3D4 element.
 *
 * These coordinates define the reference configuration of the element in
 * its local coordinate system.
 *
 * Node coordinates are:
 * - Node 1: (1, 0, 0)
 * - Node 2: (0, 1, 0)
 * - Node 3: (0, 0, 0)
 * - Node 4: (0, 0, 1)
 *
 * @return StaticMatrix<4, 3> The local coordinates of the nodes.
 */
StaticMatrix<4, 3> C3D4::node_coords_local() {
    StaticMatrix<4, 3> res {};
    res.setZero();

    // Define local coordinates for each node
    res(0, 0) = 1;  // Node 1: (1, 0, 0)
    res(1, 1) = 1;  // Node 2: (0, 1, 0)
    res(2, 2) = 0;  // Node 3: (0, 0, 0)
    res(3, 2) = 1;  // Node 4: (0, 0, 1)

    return res;
}

}  // namespace model
}  // namespace fem
