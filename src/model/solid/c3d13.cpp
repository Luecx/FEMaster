/**
 * @file c3d13.cpp
 * @brief Implements the thirteen-node quadratic pyramid solid geometry.
 *
 * Rational interpolation maps the natural pyramid, including its collapsed
 * apex. Boundary surface extraction is currently unavailable.
 * The topology supplies natural interpolation, node ordering and quadrature
 * to SolidElement; material response and assembly belong to the common base.
 *
 * Natural-space polynomial recovery is cached in extrapolation_matrix();
 * the common solid formulation applies it to integration-point output.
 *
 * @see SolidElement
 */

#include "c3d13.h"

#include "../../math/extrapolate.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface3.h"

namespace fem {
namespace model {

C3D13::C3D13(ID p_elem_id, const std::array<ID, 13>& p_node_ids)
    : SolidElement(p_elem_id, p_node_ids) {}

/**
 * Returns the static natural-domain volume integration rule.
 *
 * Mass, volume, distributed fields and thermal operators use this rule. The
 * stiffness rule may be selected separately for constitutive integration.
 *
 * @return Shared immutable topology quadrature.
 */
const math::quadrature::Quadrature& C3D13::integration_scheme() const {
    const static math::quadrature::Quadrature quad {math::quadrature::DOMAIN_ISO_PYRAMID, math::quadrature::ORDER_QUARTIC};
    return quad;
}
/**
 * Returns no boundary surface for the quadratic pyramid.
 *
 * Surface extraction is not implemented for this topology. Volume mechanics
 * remain available, but callers cannot obtain a boundary face through this API.
 *
 * @param surface_id Requested face identifier.
 * @return nullptr for every face identifier.
 */
SurfacePtr C3D13::surface(ID surface_id) {
    // TODO: implement this
    (void) surface_id;
    return nullptr;
}

/**
 * Evaluates the thirteen-node quadratic pyramid shape functions in natural coordinates.
 *
 * Rational interpolation maps the natural pyramid, including its collapsed
 * apex. Boundary surface extraction is currently unavailable.
 * Values follow the fixed connectivity ordering used by nodal interpolation.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return One interpolation weight per element node.
 */
StaticMatrix<13, 1> C3D13::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<13, 1> res {};

    // Define shape functions as provided
    res(0, 0)  = 0.25 * (r + s - 1) * ((1 + r) * (1 + s) - t + (r * s * t) / (1 - t));
    res(1, 0)  = 0.25 * (-r + s - 1) * ((1 - r) * (1 + s) - t - (r * s * t) / (1 - t));
    res(2, 0)  = 0.25 * (-r - s - 1) * ((1 - r) * (1 - s) - t + (r * s * t) / (1 - t));
    res(3, 0)  = 0.25 * (r - s - 1) * ((1 + r) * (1 - s) - t - (r * s * t) / (1 - t));
    res(4, 0)  = t * (2 * t - 1);
    res(5, 0)  = ((1 + r - t) * (1 - r - t) * (1 + s - t)) / (2 * (1 - t));
    res(6, 0)  = ((1 + s - t) * (1 - s - t) * (1 - r - t)) / (2 * (1 - t));
    res(7, 0)  = ((1 + r - t) * (1 - r - t) * (1 - s - t)) / (2 * (1 - t));
    res(8, 0)  = ((1 + s - t) * (1 - s - t) * (1 + r - t)) / (2 * (1 - t));
    res(9, 0)  = t * (1 + r - t) * (1 + s - t) / (1 - t);
    res(10, 0) = t * (1 - r - t) * (1 + s - t) / (1 - t);
    res(11, 0) = t * (1 - r - t) * (1 - s - t) / (1 - t);
    res(12, 0) = t * (1 + r - t) * (1 - s - t) / (1 - t);

    return res;
}

/**
 * Evaluates natural derivatives of the thirteen-node quadratic pyramid interpolation.
 *
 * Rows follow connectivity; columns contain dN/dr, dN/ds and dN/dt. The
 * common solid Jacobian transforms these derivatives into global gradients.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Natural shape-function derivative matrix.
 */
StaticMatrix<13, 3> C3D13::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<13, 3> der {};

    // Define shape function derivatives as provided
    der(0, 0) = 0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (r + 1) * (s + 1)
                + (0.25 * r + 0.25 * s - 0.25) * (s * t / (1 - t) + s + 1);
    der(0, 1) = 0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (r + 1) * (s + 1)
                + (0.25 * r + 0.25 * s - 0.25) * (r * t / (1 - t) + r + 1);
    der(0, 2) = (0.25 * r + 0.25 * s - 0.25) * (r * s * t / (1 - t) / (1 - t) + r * s / (1 - t) - 1);

    der(1, 0) = 0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - r) * (s + 1)
                + (-0.25 * r + 0.25 * s - 0.25) * (-s * t / (1 - t) - s - 1);
    der(1, 1) = -0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (1 - r) * (s + 1)
                + (-0.25 * r + 0.25 * s - 0.25) * (-r * t / (1 - t) - r + 1);
    der(1, 2) = (-0.25 * r + 0.25 * s - 0.25) * (-r * s * t / (1 - t) / (1 - t) - r * s / (1 - t) - 1);

    der(2, 0) = -0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - r) * (1 - s)
                + (-0.25 * r - 0.25 * s - 0.25) * (s * t / (1 - t) + s - 1);
    der(2, 1) = -0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - r) * (1 - s)
                + (-0.25 * r - 0.25 * s - 0.25) * (r * t / (1 - t) + r - 1);
    der(2, 2) = (-0.25 * r - 0.25 * s - 0.25) * (r * s * t / (1 - t) / (1 - t) + r * s / (1 - t) - 1);

    der(3, 0) = -0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (1 - s) * (r + 1)
                + (0.25 * r - 0.25 * s - 0.25) * (-s * t / (1 - t) - s + 1);
    der(3, 1) = 0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - s) * (r + 1)
                + (0.25 * r - 0.25 * s - 0.25) * (-r * t / (1 - t) - r - 1);
    der(3, 2) = (0.25 * r - 0.25 * s - 0.25) * (-r * s * t / (1 - t) / (1 - t) - r * s / (1 - t) - 1);

    der(4, 0) = 0;
    der(4, 1) = 0;
    der(4, 2) = 4 * t - 1;

    der(5, 0) = (-r - t + 1) * (s - t + 1) / (2 - 2 * t) - (r - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(5, 1) = (-r - t + 1) * (r - t + 1) / (2 - 2 * t);
    der(5, 2) = -(-r - t + 1) * (r - t + 1) / (2 - 2 * t) - (-r - t + 1) * (s - t + 1) / (2 - 2 * t)
                - (r - t + 1) * (s - t + 1) / (2 - 2 * t)
                + 2 * (-r - t + 1) * (r - t + 1) * (s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(6, 0) = -(-s - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(6, 1) = (-r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (-r - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(6, 2) = -(-r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (-r - t + 1) * (s - t + 1) / (2 - 2 * t)
                - (-s - t + 1) * (s - t + 1) / (2 - 2 * t)
                + 2 * (-r - t + 1) * (-s - t + 1) * (s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(7, 0) = (-r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (r - t + 1) * (-s - t + 1) / (2 - 2 * t);
    der(7, 1) = -(-r - t + 1) * (r - t + 1) / (2 - 2 * t);
    der(7, 2) = -(-r - t + 1) * (r - t + 1) / (2 - 2 * t) - (-r - t + 1) * (-s - t + 1) / (2 - 2 * t)
                - (r - t + 1) * (-s - t + 1) / (2 - 2 * t)
                + 2 * (-r - t + 1) * (r - t + 1) * (-s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(8, 0) = (-s - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(8, 1) = (r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (r - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(8, 2) = -(r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (r - t + 1) * (s - t + 1) / (2 - 2 * t)
                - (-s - t + 1) * (s - t + 1) / (2 - 2 * t)
                + 2 * (r - t + 1) * (-s - t + 1) * (s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(9, 0) = t * (s - t + 1) / (1 - t);
    der(9, 1) = t * (r - t + 1) / (1 - t);
    der(9, 2) = -t * (r - t + 1) / (1 - t) - t * (s - t + 1) / (1 - t)
                + t * (r - t + 1) * (s - t + 1) / ((1 - t) * (1 - t)) + (r - t + 1) * (s - t + 1) / (1 - t);

    der(10, 0) = -t * (s - t + 1) / (1 - t);
    der(10, 1) = t * (-r - t + 1) / (1 - t);
    der(10, 2) = -t * (-r - t + 1) / (1 - t) - t * (s - t + 1) / (1 - t)
                 + t * (-r - t + 1) * (s - t + 1) / ((1 - t) * (1 - t)) + (-r - t + 1) * (s - t + 1) / (1 - t);

    der(11, 0) = -t * (-s - t + 1) / (1 - t);
    der(11, 1) = -t * (-r - t + 1) / (1 - t);
    der(11, 2) = -t * (-r - t + 1) / (1 - t) - t * (-s - t + 1) / (1 - t)
                 + t * (-r - t + 1) * (-s - t + 1) / ((1 - t) * (1 - t)) + (-r - t + 1) * (-s - t + 1) / (1 - t);

    der(12, 0) = t * (-s - t + 1) / (1 - t);
    der(12, 1) = -t * (r - t + 1) / (1 - t);
    der(12, 2) = -t * (r - t + 1) / (1 - t) - t * (-s - t + 1) / (1 - t)
                 + t * (r - t + 1) * (-s - t + 1) / ((1 - t) * (1 - t)) + (r - t + 1) * (-s - t + 1) / (1 - t);

    return der;
}

/**
 * Returns natural node positions in the fixed interpolation order.
 *
 * The same ordering is used by shape functions, global connectivity and the
 * constant integration-point extrapolation operator.
 *
 * @return One natural r, s, t coordinate row per node.
 */
StaticMatrix<13, 3> C3D13::node_coords_local() {
    StaticMatrix<13, 3> res {};
    res <<   1  ,  1  , 0,
            -1  ,  1  , 0,
            -1  , -1  , 0,
             1  , -1  , 0,
             0  ,  0  , 1,
             0  ,  1  , 0,
            -1  ,  0  , 0,
             0  , -1  , 0,
             1  ,  0  , 0,
             0.5,  0.5, 0.5,
            -0.5,  0.5, 0.5,
            -0.5, -0.5, 0.5,
             0.5, -0.5, 0.5;

    return res;
}

/**
 * @brief Returns the cached natural-space integration-point-to-node recovery map.
 *
 * Rows of the operator correspond to node_coords_local() order; columns follow
 * stress_strain_ip_rst() and therefore the constitutive quadrature/state order.
 * For one scalar component q, nodal recovery is q_node = E * q_ip. Vector and
 * tensor components use the same operator independently in their existing basis.
 *
 * The ten monomials span a complete quadratic field in the natural coordinates.
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
const RowMatrix& C3D13::extrapolation_matrix() {
    // Fit the topology-selected polynomial basis at constitutive points and evaluate
    // it at the natural node coordinates. The static operator is reused by all instances.
    static const RowMatrix matrix = math::extrapolate(
        this->stress_strain_ip_rst(), this->node_coords_local(),
        {math::ExtrapolationBasis::F1,
         math::ExtrapolationBasis::FR,
         math::ExtrapolationBasis::FS,
         math::ExtrapolationBasis::FT,
         math::ExtrapolationBasis::FRR,
         math::ExtrapolationBasis::FSS,
         math::ExtrapolationBasis::FTT,
         math::ExtrapolationBasis::FRS,
         math::ExtrapolationBasis::FRT,
         math::ExtrapolationBasis::FST});
    // Return the cached node-by-integration-point map without changing element state.
    return matrix;
}

}    // namespace model
}    // namespace fem
