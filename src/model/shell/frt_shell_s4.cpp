/**
 * @file frt_shell_s4.cpp
 * @brief Implements the four-node finite-rotation MITC4 shell element.
 *
 * MITC4 interpolation is performed in natural covariant components before the
 * common shell kernel transforms the final generalized strains into the
 * pointwise orthonormal reference basis. Values and B rows use identical tying
 * coefficients; geometric-resultant weights use the exact transpose.
 *
 * @see FRTShellS4
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell_s4.h"

namespace fem::model {

/**
 * Constructs the four-node MITC4 finite-rotation shell.
 *
 * @param id Element identifier.
 * @param nodes Four nodal identifiers in positive quadrilateral ordering.
 */
FRTShellS4::FRTShellS4(ID id, const std::array<ID, 4>& nodes)
    : FRTShell<4>       (id, nodes),
      integration_scheme_(math::quadrature::Domain::DOMAIN_ISO_QUAD,
                          math::quadrature::Order::ORDER_CUBIC) {}

/**
 * Returns the public element type identifier.
 *
 * @return String identifier `MITC4FRT`.
 */
std::string FRTShellS4::type_name() const {
    return "MITC4FRT";
}

/**
 * Returns an oriented four-node surface representation.
 *
 * @param surface_id Requested positive or negative shell side.
 * @return Oriented bilinear surface object.
 */
std::shared_ptr<SurfaceInterface> FRTShellS4::surface(int surface_id) {
    return std::make_shared<Surface4>(
        surface_id == 1
            ? std::array<ID, 4>{
                  this->nodes()[0], this->nodes()[1],
                  this->nodes()[2], this->nodes()[3]
              }
            : std::array<ID, 4>{
                  this->nodes()[3], this->nodes()[2],
                  this->nodes()[1], this->nodes()[0]
              }
    );
}

/**
 * Returns the numerical area quadrature used by the element.
 *
 * @return Cubic quadrilateral quadrature rule.
 */
const math::quadrature::Quadrature& FRTShellS4::integration_scheme() const {
    return integration_scheme_;
}

/**
 * Returns all numerical integration-point coordinates for result storage.
 *
 * @return Matrix of `(r,s,t)` integration-point coordinates.
 */
RowMatrix FRTShellS4::stress_strain_ip_rst() {
    RowMatrix rst(integration_scheme_.count(), 3);
    rst.setZero();

    // Preserve the native quadrature-point ordering used by element assembly
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);
        rst(ip, 0) = point.r;
        rst(ip, 1) = point.s;
    }

    return rst;
}

/**
 * Returns the four natural corner coordinates for nodal result recovery.
 *
 * @return Matrix of `(r,s,t)` nodal coordinates.
 */
RowMatrix FRTShellS4::stress_strain_nodal_rst() {
    RowMatrix rst(4, 3);
    rst << Precision(-1), Precision(-1), Precision(0),
           Precision( 1), Precision(-1), Precision(0),
           Precision( 1), Precision( 1), Precision(0),
           Precision(-1), Precision( 1), Precision(0);
    return rst;
}

/**
 * Evaluates the four bilinear quadrilateral shape functions.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Four bilinear shape-function values.
 */
FRTShellS4::VecN FRTShellS4::shape_function(Precision r, Precision s) const {
    VecN shape;
    shape(0) = Precision(0.25) * (Precision(1) - r) * (Precision(1) - s);
    shape(1) = Precision(0.25) * (Precision(1) + r) * (Precision(1) - s);
    shape(2) = Precision(0.25) * (Precision(1) + r) * (Precision(1) + s);
    shape(3) = Precision(0.25) * (Precision(1) - r) * (Precision(1) + s);
    return shape;
}

/**
 * Evaluates the first natural derivatives of the bilinear shape functions.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Matrix whose rows contain `[dN_i/dr,dN_i/ds]`.
 */
FRTShellS4::MatN2 FRTShellS4::shape_derivative(Precision r, Precision s) const {
    MatN2 derivative;

    derivative(0, 0) = -Precision(0.25) * (Precision(1) - s);
    derivative(0, 1) = -Precision(0.25) * (Precision(1) - r);
    derivative(1, 0) =  Precision(0.25) * (Precision(1) - s);
    derivative(1, 1) = -Precision(0.25) * (Precision(1) + r);
    derivative(2, 0) =  Precision(0.25) * (Precision(1) + s);
    derivative(2, 1) =  Precision(0.25) * (Precision(1) + r);
    derivative(3, 0) = -Precision(0.25) * (Precision(1) + s);
    derivative(3, 1) =  Precision(0.25) * (Precision(1) - r);

    return derivative;
}

/**
 * Returns the natural coordinates of the four quadrilateral corners.
 *
 * @return Matrix containing the node coordinates on `[-1,1]^2`.
 */
FRTShellS4::MatN2 FRTShellS4::node_coords_natural() const {
    MatN2 coordinates;
    coordinates << Precision(-1), Precision(-1),
                   Precision( 1), Precision(-1),
                   Precision( 1), Precision( 1),
                   Precision(-1), Precision( 1);
    return coordinates;
}

/**
 * Returns the four classical MITC4 edge-midpoint tying positions.
 *
 * The ordering is bottom, top, left and right.
 *
 * @return Ordered natural tying-point coordinates.
 */
std::vector<Vec2> FRTShellS4::tying_point_coordinates() const {
    return {
        Vec2(Precision( 0), Precision(-1)),
        Vec2(Precision( 0), Precision( 1)),
        Vec2(Precision(-1), Precision( 0)),
        Vec2(Precision( 1), Precision( 0))
    };
}

/**
 * Applies the classical MITC4 transverse-shear interpolation.
 *
 * `gamma_r3` is interpolated between the bottom and top tying points with `s`.
 * `gamma_s3` is interpolated between the left and right tying points with `r`.
 * The same interpolation is applied to the corresponding B rows.
 *
 * @param data Active tying-point strain and B data.
 * @param point Target natural shell point.
 * @param strain_nat Compatible natural strain vector to modify.
 * @param B_nat Optional compatible natural B matrix to modify.
 */
void FRTShellS4::apply_mitc_natural(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat
) const {
    constexpr Index gamma_r3 = 6;
    constexpr Index gamma_s3 = 7;

    const Precision bottom = Precision(0.5) * (Precision(1) - point.s);
    const Precision top    = Precision(0.5) * (Precision(1) + point.s);
    const Precision left   = Precision(0.5) * (Precision(1) - point.r);
    const Precision right  = Precision(0.5) * (Precision(1) + point.r);

    strain_nat(gamma_r3) = bottom * data.tying_strain_nat[0](gamma_r3)
                         + top    * data.tying_strain_nat[1](gamma_r3);
    strain_nat(gamma_s3) = left   * data.tying_strain_nat[2](gamma_s3)
                         + right  * data.tying_strain_nat[3](gamma_s3);

    if (B_nat) {
        // Apply identical edge interpolation to the first strain derivatives
        B_nat->row(gamma_r3) = bottom * data.tying_B_nat[0].row(gamma_r3)
                             + top    * data.tying_B_nat[1].row(gamma_r3);
        B_nat->row(gamma_s3) = left   * data.tying_B_nat[2].row(gamma_s3)
                             + right  * data.tying_B_nat[3].row(gamma_s3);
    }
}

/**
 * Applies the transpose of the classical MITC4 shear interpolation.
 *
 * @param point Integration point defining the four edge coefficients.
 * @param assumed_weights Generalized natural weights after local-basis pull-back.
 * @param compatible_weights Compatible integration-point weights to increment.
 * @param tying_weights Zeroed compatible tying-point weights to increment.
 */
void FRTShellS4::pull_back_mitc_resultants(
    const ReferencePoint& point,
    const Vec8&           assumed_weights,
    Vec8&                 compatible_weights,
    Span<Vec8>            tying_weights
) const {
    constexpr Index gamma_r3 = 6;
    constexpr Index gamma_s3 = 7;

    const Precision bottom = Precision(0.5) * (Precision(1) - point.s);
    const Precision top    = Precision(0.5) * (Precision(1) + point.s);
    const Precision left   = Precision(0.5) * (Precision(1) - point.r);
    const Precision right  = Precision(0.5) * (Precision(1) + point.r);

    // Membrane and curvature components remain compatible at the integration point
    compatible_weights.template segment<6>(0) +=
        assumed_weights.template segment<6>(0);

    // Distribute the two assumed shear weights to their edge tying values
    tying_weights[0](gamma_r3) += bottom * assumed_weights(gamma_r3);
    tying_weights[1](gamma_r3) += top    * assumed_weights(gamma_r3);
    tying_weights[2](gamma_s3) += left   * assumed_weights(gamma_s3);
    tying_weights[3](gamma_s3) += right  * assumed_weights(gamma_s3);
}

} // namespace fem::model
