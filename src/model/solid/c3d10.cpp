/**
 * @file c3d10.cpp
 * @brief Implements the ten-node quadratic tetrahedral solid geometry.
 *
 * Natural-coordinate shape functions, derivatives, node positions and
 * boundary surfaces define the topology used by SolidElement. The element
 * integration rules supply volume quadrature; material response and
 * mechanical assembly remain in the common solid evaluate() path.
 *
 * Natural-space polynomial recovery is cached in extrapolation_matrix();
 * the common solid formulation applies it to integration-point output.
 *
 * @see SolidElement
 */

#include "c3d10.h"

#include "../../math/extrapolate.h"
#include "../geometry/surface/surface6.h"

namespace fem {
namespace model {

C3D10::C3D10(ID pElemId, const std::array<ID, 10>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

/**
 * Evaluates the ten-node quadratic tetrahedral shape functions in natural coordinates.
 *
 * Quadratic tetrahedral interpolation uses four vertex and six edge nodes.
 * Volume and material quadrature are selected separately.
 * Values follow the fixed connectivity ordering used by nodal interpolation.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return One interpolation weight per element node.
 */
StaticMatrix<10, 1> C3D10::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<10, 1> res {};

    // Vertex nodes
    res(0) = (1 - r - s - t) * (2*(1 - r - s - t) - 1);
    res(1) = r * (2*r - 1);
    res(2) = s * (2*s - 1);
    res(3) = t * (2*t - 1);

    // Mid-edge nodes
    res(4) = 4 * r * (1 - r - s - t);
    res(5) = 4 * r * s;
    res(6) = 4 * s * (1 - r - s - t);
    res(7) = 4 * t * (1 - r - s - t);
    res(8) = 4 * r * t;
    res(9) = 4 * s * t;

    return res;
}

/**
 * Evaluates natural derivatives of the ten-node quadratic tetrahedral interpolation.
 *
 * Rows follow connectivity; columns contain dN/dr, dN/ds and dN/dt. The
 * common solid Jacobian transforms these derivatives into global gradients.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Natural shape-function derivative matrix.
 */
StaticMatrix<10, 3> C3D10::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<10, 3> res {};
    res.setZero();

    Precision u = 1 - r - s - t;

    // Derivatives with respect to r
    res(0, 0) = -4*u + 1;
    res(1, 0) = 4*r - 1;
    res(2, 0) = 0;
    res(3, 0) = 0;
    res(4, 0) = 4*u - 4*r;
    res(5, 0) = 4*s;
    res(6, 0) = -4*s;
    res(7, 0) = -4*t;
    res(8, 0) = 4*t;
    res(9, 0) = 0;

    // Derivatives with respect to s
    res(0, 1) = -4*u + 1;
    res(1, 1) = 0;
    res(2, 1) = 4*s - 1;
    res(3, 1) = 0;
    res(4, 1) = -4*r;
    res(5, 1) = 4*r;
    res(6, 1) = 4*u - 4*s;
    res(7, 1) = -4*t;
    res(8, 1) = 0;
    res(9, 1) = 4*t;

    // Derivatives with respect to t
    res(0, 2) = -4*u + 1;
    res(1, 2) = 0;
    res(2, 2) = 0;
    res(3, 2) = 4*t - 1;
    res(4, 2) = -4*r;
    res(5, 2) = 0;
    res(6, 2) = -4*s;
    res(7, 2) = 4*u - 4*t;
    res(8, 2) = 4*r;
    res(9, 2) = 4*s;

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
StaticMatrix<10, 3> C3D10::node_coords_local() {
    StaticMatrix<10, 3> res {};
    res.setZero();

    // Vertex nodes
    res(0, 0) = 0;
    res(0, 1) = 0;
    res(0, 2) = 0;
    res(1, 0) = 1;
    res(1, 1) = 0;
    res(1, 2) = 0;
    res(2, 0) = 0;
    res(2, 1) = 1;
    res(2, 2) = 0;
    res(3, 0) = 0;
    res(3, 1) = 0;
    res(3, 2) = 1;

    // Mid-edge nodes
    res(4, 0) = 0.5;
    res(4, 1) = 0.0;
    res(4, 2) = 0;
    res(5, 0) = 0.5;
    res(5, 1) = 0.5;
    res(5, 2) = 0.0;
    res(6, 0) = 0.0;
    res(6, 1) = 0.5;
    res(6, 2) = 0.0;
    res(7, 0) = 0.0;
    res(7, 1) = 0.0;
    res(7, 2) = 0.5;
    res(8, 0) = 0.5;
    res(8, 1) = 0.0;
    res(8, 2) = 0.5;
    res(9, 0) = 0.0;
    res(9, 1) = 0.5;
    res(9, 2) = 0.5;

    return res;
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
SurfacePtr C3D10::surface(ID surface_id) {
    switch (surface_id) {
        case 1:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[0], node_ids[1], node_ids[2], node_ids[4], node_ids[5], node_ids[6]});
        case 2:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[0], node_ids[3], node_ids[1], node_ids[7], node_ids[8], node_ids[4]});
        case 3:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[1], node_ids[3], node_ids[2], node_ids[8], node_ids[9], node_ids[5]});
        case 4:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[2], node_ids[3], node_ids[0], node_ids[9], node_ids[7], node_ids[6]});
        default: return nullptr;    // Invalid surface ID
    }
}
/**
 * Returns the static natural-domain volume integration rule.
 *
 * Mass, volume, distributed fields and thermal operators use this rule. The
 * stiffness rule may be selected separately for constitutive integration.
 *
 * @return Shared immutable topology quadrature.
 */
const math::quadrature::Quadrature& C3D10::integration_scheme() const {
    const static math::quadrature::Quadrature quad {math::quadrature::DOMAIN_ISO_TET, math::quadrature::ORDER_QUARTIC};
    return quad;
}
/**
 * Returns the static quadrature used by constitutive integration.
 *
 * The point ordering defines material-state rows and integration-point output.
 * This rule may contain fewer points than the topology volume quadrature.
 *
 * @return Shared immutable material quadrature.
 */
const math::quadrature::Quadrature& C3D10::integration_scheme_stiffness() const {
    const static math::quadrature::Quadrature quad {math::quadrature::DOMAIN_ISO_TET, math::quadrature::ORDER_QUADRATIC};
    return quad;
}

/**
 * @brief Returns the cached natural-space integration-point-to-node recovery map.
 *
 * Rows of the operator correspond to node_coords_local() order; columns follow
 * stress_strain_ip_rst() and therefore the constitutive quadrature/state order.
 * For one scalar component q, nodal recovery is q_node = E * q_ip. Vector and
 * tensor components use the same operator independently in their existing basis.
 *
 * The basis (1, r, s, t) reconstructs an affine field on the natural tetrahedron.
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
const RowMatrix& C3D10::extrapolation_matrix() {
    // Fit the topology-selected polynomial basis at constitutive points and evaluate
    // it at the natural node coordinates. The static operator is reused by all instances.
    static const RowMatrix matrix = math::extrapolate(
        this->stress_strain_ip_rst(), this->node_coords_local(),
        {math::ExtrapolationBasis::F1,
         math::ExtrapolationBasis::FR,
         math::ExtrapolationBasis::FS,
         math::ExtrapolationBasis::FT});
    // Return the cached node-by-integration-point map without changing element state.
    return matrix;
}

}    // namespace model
}    // namespace fem
