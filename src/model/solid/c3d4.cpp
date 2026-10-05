/**
 * @file c3d4.cpp
 * @brief Implements the four-node linear tetrahedral solid geometry.
 *
 * The natural tetrahedron uses constant shape gradients and one material integration point.
 * The topology supplies natural interpolation, node ordering and quadrature
 * to SolidElement; material response and assembly belong to the common base.
 *
 * Natural-space polynomial recovery is cached in extrapolation_matrix();
 * the common solid formulation applies it to integration-point output.
 *
 * @see SolidElement
 */

#include "c3d4.h"

#include "../../math/extrapolate.h"
#include "../geometry/surface/surface3.h"

namespace fem {
namespace model {

/**
 * @brief Constructor that initializes the C3D4 element with an element ID
 * and the corresponding node IDs.
 *
 * @param pElemId The unique ID of the element.
 * @param pNodeIds Array containing IDs of the 4 nodes that define the element.
 */
C3D4::C3D4(ID pElemId, const std::array<ID, 4>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

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
    (void) r;
    (void) s;
    (void) t;

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
    res(0, 0) = 1;    // Node 1: (1, 0, 0)
    res(1, 1) = 1;    // Node 2: (0, 1, 0)
    res(2, 2) = 0;    // Node 3: (0, 0, 0)
    res(3, 2) = 1;    // Node 4: (0, 0, 1)

    return res;
}

/**
 * Returns the static natural-domain volume integration rule.
 *
 * Mass, volume, distributed fields and thermal operators use this rule. The
 * stiffness rule may be selected separately for constitutive integration.
 *
 * @return Shared immutable topology quadrature.
 */
const math::quadrature::Quadrature& C3D4::integration_scheme() const {
    const static math::quadrature::Quadrature quad {math::quadrature::DOMAIN_ISO_TET, math::quadrature::ORDER_LINEAR};
    return quad;
}

/**
 * Extracts a boundary face using the topology-specific face numbering.
 *
 * The surface shares global node identifiers with the solid and provides its
 * own interpolation and geometric integration for loads and constraints.
 *
 * @param surface_id One-based face identifier.
 * @return Requested face, or nullptr when the identifier is invalid.
 */
SurfacePtr C3D4::surface(ID surface_id) {
    switch (surface_id) {
        case 1: return std::make_shared<Surface3>(std::array<ID, 3> {node_ids[0], node_ids[1], node_ids[2]});
        case 2: return std::make_shared<Surface3>(std::array<ID, 3> {node_ids[0], node_ids[3], node_ids[1]});
        case 3: return std::make_shared<Surface3>(std::array<ID, 3> {node_ids[1], node_ids[3], node_ids[2]});
        case 4: return std::make_shared<Surface3>(std::array<ID, 3> {node_ids[2], node_ids[3], node_ids[0]});
        default: return nullptr;    // Invalid surface ID
    }
}

/**
 * @brief Returns the cached natural-space integration-point-to-node recovery map.
 *
 * Rows of the operator correspond to node_coords_local() order; columns follow
 * stress_strain_ip_rst() and therefore the constitutive quadrature/state order.
 * For one scalar component q, nodal recovery is q_node = E * q_ip. Vector and
 * tensor components use the same operator independently in their existing basis.
 *
 * A constant basis replicates the single integration-point value at every node.
 * math::extrapolate() fits coefficients through the source normal equations and
 * evaluates the polynomial at the natural node locations. Point-count and solver
 * validation are performed by that helper; no fallback is added here.
 *
 * Function-local static initialization builds the operator once from the fixed
 * natural topology and constitutive rule. No physical geometry, material history
 * or element runtime state is stored or changed by the recovery map.
 *
 * @return Shared immutable node-by-integration-point matrix in reference space.
 */
const RowMatrix& C3D4::extrapolation_matrix() {
    // Fit the topology-selected polynomial basis at constitutive points and evaluate
    // it at the natural node coordinates. The static operator is reused by all instances.
    static const RowMatrix matrix = math::extrapolate(
        this->stress_strain_ip_rst(), this->node_coords_local(),
        {math::ExtrapolationBasis::F1});
    // Return the cached node-by-integration-point map without changing element state.
    return matrix;
}

}  // namespace model
}  // namespace fem
