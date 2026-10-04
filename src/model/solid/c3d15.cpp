/**
 * @file c3d15.cpp
 * @brief Implements the fifteen-node quadratic wedge solid geometry.
 *
 * Quadratic interpolation uses six vertex nodes and nine edge nodes on the
 * natural triangular prism.
 * The topology supplies natural interpolation, node ordering and quadrature
 * to SolidElement; material response and assembly belong to the common base.
 *
 * @see SolidElement
 */

#include "c3d15.h"

#include "../geometry/surface/surface8.h"
#include "../geometry/surface/surface6.h"

namespace fem {
namespace model {

C3D15::C3D15(ID p_elem_id, const std::array<ID, 15>& p_node_ids)
    : SolidElement(p_elem_id, p_node_ids) {}
/**
 * Returns the static natural-domain volume integration rule.
 *
 * Mass, volume, distributed fields and thermal operators use this rule. The
 * stiffness rule may be selected separately for constitutive integration.
 *
 * @return Shared immutable topology quadrature.
 */
const math::quadrature::Quadrature& C3D15::integration_scheme() const {
    const static math::quadrature::Quadrature quad {math::quadrature::DOMAIN_ISO_WEDGE, math::quadrature::ORDER_SUPER_QUADRATIC};
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
SurfacePtr C3D15::surface(ID surface_id) {
    // Quadratic wedge faces use six-node triangles and eight-node quadrilaterals
    switch (surface_id) {
        case 1:
            return std::make_shared<Surface6>(std::array<ID, 6> {node_ids[0],
                                                                 node_ids[1],
                                                                 node_ids[2],
                                                                 node_ids[6],
                                                                 node_ids[7],
                                                                 node_ids[8]});    // Face 1: Triangle
        case 2:
            return std::make_shared<Surface6>(std::array<ID, 6> {node_ids[3],
                                                                 node_ids[4],
                                                                 node_ids[5],
                                                                 node_ids[10],
                                                                 node_ids[11],
                                                                 node_ids[12]});    // Face 2: Triangle
        case 3:
            return std::make_shared<Surface8>(std::array<ID, 8> {node_ids[0],
                                                                 node_ids[1],
                                                                 node_ids[4],
                                                                 node_ids[3],
                                                                 node_ids[6],
                                                                 node_ids[13],
                                                                 node_ids[9],
                                                                 node_ids[12]});    // Face 3: Quadrilateral
        case 4:
            return std::make_shared<Surface8>(std::array<ID, 8> {node_ids[1],
                                                                 node_ids[2],
                                                                 node_ids[5],
                                                                 node_ids[4],
                                                                 node_ids[7],
                                                                 node_ids[14],
                                                                 node_ids[10],
                                                                 node_ids[13]});    // Face 4: Quadrilateral
        case 5:
            return std::make_shared<Surface8>(std::array<ID, 8> {node_ids[2],
                                                                 node_ids[0],
                                                                 node_ids[3],
                                                                 node_ids[5],
                                                                 node_ids[8],
                                                                 node_ids[12],
                                                                 node_ids[11],
                                                                 node_ids[14]});    // Face 5: Quadrilateral
        default: return nullptr;                                                    // Invalid surface ID
    }
}

/**
 * Evaluates the fifteen-node quadratic wedge shape functions in natural coordinates.
 *
 * Quadratic interpolation uses six vertex nodes and nine edge nodes on the
 * natural triangular prism.
 * Values follow the fixed connectivity ordering used by nodal interpolation.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return One interpolation weight per element node.
 */
StaticMatrix<15, 1> C3D15::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<15, 1> res {};

    Precision           g = r;
    Precision           h = s;
    r                     = t;

    Precision one_mg_mh = 1 - g - h;
    Precision one_mr    = 1 - r;
    Precision one_pr    = 1 + r;
    Precision one_mr2   = 1 - r * r;
    Precision two_g_m1  = 2 * g - 1;
    Precision two_h_m1  = 2 * h - 1;

    // shape function values
    res(0 ,0) = 0.5 * (one_mg_mh * (2 * one_mg_mh - 1) * one_mr - one_mg_mh * one_mr2);
    res(1 ,0) = 0.5 * (g * two_g_m1 * one_mr - g * one_mr2);
    res(2 ,0) = 0.5 * (h * two_h_m1 * one_mr - h * one_mr2);
    res(3 ,0) = 0.5 * (one_mg_mh * (2 * one_mg_mh - 1) * one_pr - one_mg_mh * one_mr2);
    res(4 ,0) = 0.5 * (g * two_g_m1 * one_pr - g * one_mr2);
    res(5 ,0) = 0.5 * (h * two_h_m1 * one_pr - h * one_mr2);
    res(6 ,0) = 2 * one_mg_mh * g * one_mr;
    res(7 ,0) = 2 * g * h * one_mr;
    res(8 ,0) = 2 * h * one_mg_mh * one_mr;
    res(9 ,0) = 2 * one_mg_mh * g * one_pr;
    res(10,0) = 2 * g * h * one_pr;
    res(11,0) = 2 * h * one_mg_mh * one_pr;
    res(12,0) = one_mg_mh * one_mr2;
    res(13,0) = g * one_mr2;
    res(14,0) = h * one_mr2;

    return res;
}

/**
 * Evaluates natural derivatives of the fifteen-node quadratic wedge interpolation.
 *
 * Rows follow connectivity; columns contain dN/dr, dN/ds and dN/dt. The
 * common solid Jacobian transforms these derivatives into global gradients.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Natural shape-function derivative matrix.
 */
StaticMatrix<15, 3> C3D15::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<15, 3> res {};
    res.setZero();

    // Derivative with respect to g
    res(0, 0) = -0.5 * (-1 + t) * (-2 + 4 * r + 4 * s + t);
    res(0, 1) = -0.5 * (-1 + t) * (-2 + 4 * r + 4 * s + t);
    res(0, 2) = -0.5 * (-1 + r + s) * (-1 + 2 * r + 2 * s + 2 * t);

    res(1, 0) = -1 + r * (2 - 2 * t) + 0.5 * t + 0.5 * t * t;
    res(1, 1) = 0;
    res(1, 2) = r * (0.5 - r + t);

    res(2, 0) = 0;
    res(2, 1) = -1 + s * (2 - 2 * t) + 0.5 * t + 0.5 * t * t;
    res(2, 2) = s * (0.5 - s + t);

    res(3, 0) = 0.5 * (1 + t) * (-2 + 4 * r + 4 * s - t);
    res(3, 1) = 0.5 * (1 + t) * (-2 + 4 * r + 4 * s - t);
    res(3, 2) = 0.5 * (-1 + r + s) * (-1 + 2 * r + 2 * s - 2 * t);

    res(4, 0) = -1 - 0.5 * t + 0.5 * t * t + r * (2 + 2 * t);
    res(4, 1) = 0;
    res(4, 2) = 0.5 * r * (-1 + 2 * r + 2 * t);

    res(5, 0) = 0;
    res(5, 1) = -1 - 0.5 * t + 0.5 * t * t + s * (2 + 2 * t);
    res(5, 2) = 0.5 * s * (-1 + 2 * s + 2 * t);

    res(6, 0) = (1 - t) * (2 * (1 - r - s) - 2 * r);
    res(6, 1) = -2 * r * (1 - t);
    res(6, 2) = -2 * r * (1 - r - s);

    res(7, 0) = 2 * s * (1 - t);
    res(7, 1) = 2 * r * (1 - t);
    res(7, 2) = -2 * r * s;

    res(8, 0) = -2 * s * (1 - t);
    res(8, 1) = (1 - t) * (2 * (1 - r - s) - 2 * s);
    res(8, 2) = -2 * s * (1 - r - s);

    res(9, 0) = (1 + t) * (2 * (1 - r - s) - 2 * r);
    res(9, 1) = -2 * r * (1 + t);
    res(9, 2) = 2 * r * (1 - r - s);

    res(10, 0) = 2 * s * (1 + t);
    res(10, 1) = 2 * r * (1 + t);
    res(10, 2) = 2 * r * s;

    res(11, 0) = -2 * s * (1 + t);
    res(11, 1) = (1 + t) * (2 * (1 - r - s) - 2 * s);
    res(11, 2) = 2 * s * (1 - r - s);

    res(12, 0) = -(1 - t * t);
    res(12, 1) = -(1 - t * t);
    res(12, 2) = (1 - r - s) * (-2 * t);

    res(13, 0) = (1 - t * t);
    res(13, 1) = 0;
    res(13, 2) = -r * 2 * t;

    res(14, 0) = 0;
    res(14, 1) = (1 - t * t);
    res(14, 2) = -s * 2 * t;

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
StaticMatrix<15, 3> C3D15::node_coords_local() {
    StaticMatrix<15, 3> res {};
    res.setZero();

    // Vertex nodes
    res(0, 0) = 0;
    res(0, 1) = 0;
    res(0, 2) = -1;    // Node 3
    res(1, 0) = 1;
    res(1, 1) = 0;
    res(1, 2) = -1;    // Node 1
    res(2, 0) = 0;
    res(2, 1) = 1;
    res(2, 2) = -1;    // Node 2
    res(3, 0) = 0;
    res(3, 1) = 0;
    res(3, 2) = 1;    // Node 6
    res(4, 0) = 1;
    res(4, 1) = 0;
    res(4, 2) = 1;    // Node 4
    res(5, 0) = 0;
    res(5, 1) = 1;
    res(5, 2) = 1;    // Node 5

    // Mid-edge nodes (bottom face, then top face)
    res(6, 0)  = 0.5;
    res(6, 1)  = 0;
    res(6, 2)  = -1;    // Node 9
    res(7, 0)  = 0.5;
    res(7, 1)  = 0.5;
    res(7, 2)  = -1;    // Node 7
    res(8, 0)  = 0;
    res(8, 1)  = 0.5;
    res(8, 2)  = -1;    // Node 8
    res(9, 0)  = 0.5;
    res(9, 1)  = 0;
    res(9, 2)  = 1;    // Node 12
    res(10, 0) = 0.5;
    res(10, 1) = 0.5;
    res(10, 2) = 1;    // Node 10
    res(11, 0) = 0;
    res(11, 1) = 0.5;
    res(11, 2) = 1;    // Node 11

    // Mid-edge nodes (connecting top and bottom faces)
    res(12, 0) = 0;
    res(12, 1) = 0.0;
    res(12, 2) = 0;    // Node 15
    res(13, 0) = 1.0;
    res(13, 1) = 0.0;
    res(13, 2) = 0;    // Node 13
    res(14, 0) = 0.0;
    res(14, 1) = 1.0;
    res(14, 2) = 0;    // Node 14

    return res;
}
}    // namespace model
}    // namespace fem
