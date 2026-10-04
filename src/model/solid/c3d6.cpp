/**
 * @file c3d6.cpp
 * @brief Implements the six-node linear wedge solid geometry.
 *
 * The natural triangular prism combines linear triangular interpolation with linear interpolation
 * through its thickness.
 * The topology supplies natural interpolation, node ordering and quadrature
 * to SolidElement; material response and assembly belong to the common base.
 *
 * @see SolidElement
 */

#include "c3d6.h"

#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface3.h"

namespace fem {
namespace model {

C3D6::C3D6(ID p_elem_id, const std::array<ID, 6>& p_node_ids)
    : SolidElement(p_elem_id, p_node_ids) {}

/**
 * Returns the static natural-domain volume integration rule.
 *
 * Mass, volume, distributed fields and thermal operators use this rule. The
 * stiffness rule may be selected separately for constitutive integration.
 *
 * @return Shared immutable topology quadrature.
 */
const math::quadrature::Quadrature& C3D6::integration_scheme() const {
    const static math::quadrature::Quadrature quad {math::quadrature::DOMAIN_ISO_WEDGE, math::quadrature::ORDER_SUPER_LINEAR};
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
SurfacePtr C3D6::surface(ID surface_id) {
    // C3D6 (6-node wedge element): Triangular faces have 3 nodes, quadrilateral faces have 4 nodes
    switch (surface_id) {
        case 1:
            return std::make_shared<Surface3>(
                std::array<ID, 3> {node_ids[0], node_ids[1], node_ids[2]});    // Face 1: Triangle
        case 2:
            return std::make_shared<Surface3>(
                std::array<ID, 3> {node_ids[3], node_ids[4], node_ids[5]});    // Face 2: Triangle
        case 3:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[0], node_ids[1], node_ids[4], node_ids[3]});    // Face 3: Quadrilateral
        case 4:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[1], node_ids[2], node_ids[5], node_ids[4]});    // Face 4: Quadrilateral
        case 5:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[2], node_ids[0], node_ids[3], node_ids[5]});    // Face 5: Quadrilateral
        default: return nullptr;                                                            // Invalid surface ID
    }
}

/**
 * Evaluates the six-node linear wedge shape functions in natural coordinates.
 *
 * The natural triangular prism combines linear triangular interpolation with linear interpolation
 * through its thickness.
 * Values follow the fixed connectivity ordering used by nodal interpolation.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return One interpolation weight per element node.
 */
StaticMatrix<6, 1> C3D6::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<6, 1> res {};

    // g = r
    // h = s
    // r = t

    // Vertex nodes
    res(0) = r * (1 - t) / 2;
    res(1) = s * (1 - t) / 2;
    res(2) = (1 - r - s) * (1 - t) / 2;
    res(3) = r * (1 + t) / 2;
    res(4) = s * (1 + t) / 2;
    res(5) = (1 - r - s) * (1 + t) / 2;

    return res;
}

/**
 * Evaluates natural derivatives of the six-node linear wedge interpolation.
 *
 * Rows follow connectivity; columns contain dN/dr, dN/ds and dN/dt. The
 * common solid Jacobian transforms these derivatives into global gradients.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Natural shape-function derivative matrix.
 */
StaticMatrix<6, 3> C3D6::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<6, 3> res {};
    res.setZero();

    // Derivatives with respect to r
    res(0, 0) = (1 - t) / 2;
    res(1, 0) = 0;
    res(2, 0) = -(1 - t) / 2;
    res(3, 0) = (1 + t) / 2;
    res(4, 0) = 0;
    res(5, 0) = -(1 + t) / 2;

    // Derivatives with respect to s
    res(0, 1) = 0;
    res(1, 1) = (1 - t) / 2;
    res(2, 1) = -(1 - t) / 2;
    res(3, 1) = 0;
    res(4, 1) = (1 + t) / 2;
    res(5, 1) = -(1 + t) / 2;

    // Derivatives with respect to t
    res(0, 2) = -r / 2;
    res(1, 2) = -s / 2;
    res(2, 2) = -(1 - r - s) / 2;
    res(3, 2) = r / 2;
    res(4, 2) = s / 2;
    res(5, 2) = (1 - r - s) / 2;

    return res;
}

/**
 * Returns natural node positions in the fixed interpolation order.
 *
 * The same ordering is used by shape functions, global connectivity and the
 * constant integration-point extrapolation operator.
 *
 * @return One natural r, s, t coordinate row per node.
 */
StaticMatrix<6, 3> C3D6::node_coords_local() {
    StaticMatrix<6, 3> res {};
    res.setZero();

    // Vertex nodes
    res(0, 0) = 1;   res(0, 1) = 0;   res(0, 2) = -1;  // Node 1
    res(1, 0) = 0;   res(1, 1) = 1;   res(1, 2) = -1;  // Node 2
    res(2, 0) = 0;   res(2, 1) = 0;   res(2, 2) = -1;  // Node 3
    res(3, 0) = 1;   res(3, 1) = 0;   res(3, 2) = 1;   // Node 4
    res(4, 0) = 0;   res(4, 1) = 1;   res(4, 2) = 1;   // Node 5
    res(5, 0) = 0;   res(5, 1) = 0;   res(5, 2) = 1;   // Node 6

    return res;
}
}    // namespace model
}    // namespace fem
