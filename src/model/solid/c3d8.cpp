/**
 * @file c3d8.cpp
 * @brief Implements the eight-node trilinear hexahedral solid element geometry.
 *
 * `C3D8` provides the natural-coordinate interpolation, natural derivatives,
 * reference-node locations, face extraction and integration rule required by the
 * common solid-element mechanics implemented in `SolidElement`.
 *
 * Surface extraction follows the element's six Abaqus-style face ids. Each face
 * is returned with an outward node winding so surface normals are consistent for
 * pressure, tie and contact formulations. In particular, S1 is the `t = -1`
 * face and therefore uses the opposite winding from the parallel S2 face.
 *
 * @see C3D8
 * @see SolidElement
 * @see Surface4
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#include "c3d8.h"

#include "../geometry/surface/surface4.h"

namespace fem {
namespace model {

/**
 * Constructs one C3D8 element from its global element id and eight global node
 * ids. Kinematic and constitutive behavior is inherited from `SolidElement`.
 */
C3D8::C3D8(ID pElemId, const std::array<ID, 8>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

/**
 * Evaluates the eight trilinear C3D8 shape functions in natural coordinates.
 *
 * The local coordinates `(r, s, t)` lie in `[-1, 1]^3`. The implementation
 * constructs the signed factors `(1 +/- r)(1 +/- s)(1 +/- t)` for the element's
 * fixed node ordering and applies the common factor `1/8`.
 *
 * @return Column vector containing `N_1 ... N_8`.
 */
StaticMatrix<8, 1> C3D8::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<8, 1> res {};

    const Precision rp = r + Precision(1);
    const Precision sp = s + Precision(1);
    const Precision tp = t + Precision(1);
    const Precision rm = r - Precision(1);
    const Precision sm = s - Precision(1);
    const Precision tm = t - Precision(1);

    // Each node is represented by one sign and one factor in r, s and t.
    for (int n = 0; n < 8; ++n) {
        const Precision sign = (n == 0 || n == 2 || n == 5 || n == 7) ? Precision(-1) : Precision(1);
        const Precision r1   = ((n + 1) / 2) % 2 == 1 ? rp : rm;
        const Precision s1   = (n / 2) % 2 == 1 ? sp : sm;
        const Precision t1   = n >= 4 ? tp : tm;

        res(n, 0) = sign * r1 * s1 * t1;
    }

    res *= Precision(0.125);
    return res;
}

/**
 * Evaluates natural derivatives of the eight trilinear shape functions.
 *
 * Columns correspond to derivatives with respect to `(r, s, t)`. For each
 * derivative, the factor belonging to the differentiated coordinate is removed
 * from the trilinear product while the remaining two factors and node sign are
 * retained. The common `1/8` scale is applied after assembly.
 *
 * @return `8 x 3` matrix containing `(dN/dr, dN/ds, dN/dt)` per node.
 */
StaticMatrix<8, 3> C3D8::shape_derivative(Precision r, Precision s, Precision t) {
    const Precision rp = r + Precision(1);
    const Precision sp = s + Precision(1);
    const Precision tp = t + Precision(1);
    const Precision rm = r - Precision(1);
    const Precision sm = s - Precision(1);
    const Precision tm = t - Precision(1);

    // Cache the sign and three natural-coordinate factors for all nodes.
    StaticMatrix<8, 4> factors {};
    for (int n = 0; n < 8; ++n) {
        factors(n, 0) = (n == 0 || n == 2 || n == 5 || n == 7) ? Precision(-1) : Precision(1);
        factors(n, 1) = ((n + 1) / 2) % 2 == 1 ? rp : rm;
        factors(n, 2) = (n / 2) % 2 == 1 ? sp : sm;
        factors(n, 3) = n >= 4 ? tp : tm;
    }

    StaticMatrix<8, 3> derivative {};
    for (int n = 0; n < 8; ++n) {
        derivative(n, 0) = factors(n, 0) * factors(n, 2) * factors(n, 3);
        derivative(n, 1) = factors(n, 0) * factors(n, 1) * factors(n, 3);
        derivative(n, 2) = factors(n, 0) * factors(n, 1) * factors(n, 2);
    }

    derivative *= Precision(0.125);
    return derivative;
}

/**
 * Returns the natural coordinates of the eight C3D8 nodes.
 *
 * Node coordinates follow the same fixed ordering used by the interpolation and
 * face definitions, with each coordinate taking the value `-1` or `+1`.
 */
StaticMatrix<8, 3> C3D8::node_coords_local() {
    StaticMatrix<8, 3> res {};

    for (int n = 0; n < 8; ++n) {
        const Precision r = ((n + 1) / 2) % 2 == 0 ? Precision(-1) : Precision(1);
        const Precision s = (n / 2) % 2 == 0 ? Precision(-1) : Precision(1);
        const Precision t = n >= 4 ? Precision(1) : Precision(-1);

        res(n, 0) = r;
        res(n, 1) = s;
        res(n, 2) = t;
    }

    return res;
}

/**
 * Extracts one quadrilateral C3D8 boundary face with outward node winding.
 *
 * Surface ids follow the existing solid-face convention. The returned `Surface4`
 * reuses the element's global node ids. S1 is the lower `t = -1` face and must be
 * wound opposite to S2 so its computed normal points outward in the `-t`
 * direction. The remaining four side faces are ordered analogously.
 *
 * @param surface_id Face id in the range `1 ... 6`.
 * @return Surface representation of the requested face, or `nullptr` for an
 *         invalid id.
 */
SurfacePtr C3D8::surface(ID surface_id) {
    switch (surface_id) {
        case 1:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[0], node_ids[3], node_ids[2], node_ids[1]});
        case 2:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[4], node_ids[5], node_ids[6], node_ids[7]});
        case 3:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[0], node_ids[1], node_ids[5], node_ids[4]});
        case 4:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[1], node_ids[2], node_ids[6], node_ids[5]});
        case 5:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[2], node_ids[3], node_ids[7], node_ids[6]});
        case 6:
            return std::make_shared<Surface4>(
                std::array<ID, 4> {node_ids[3], node_ids[0], node_ids[4], node_ids[7]});
        default:
            return nullptr;
    }
}

/**
 * Returns the element's full quadratic-order hexahedral integration scheme.
 *
 * The quadrature object is static because all C3D8 elements use the same natural
 * integration rule and the rule itself carries no element-specific state.
 */
const math::quadrature::Quadrature& C3D8::integration_scheme() const {
    static const math::quadrature::Quadrature quadrature {
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };
    return quadrature;
}

} // namespace model
} // namespace fem
