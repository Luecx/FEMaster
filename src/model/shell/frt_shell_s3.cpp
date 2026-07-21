/**
 * @file frt_shell_s3.cpp
 * @brief Implements the three-node finite-rotation MITC3 shell element.
 *
 * The assumed covariant shear field follows the isotropic MITC3 interpolation
 * on the natural triangle. Identical coefficients are applied to shear values
 * and B rows, while generalized resultant weights use the exact transposed
 * interpolation during geometric-tangent assembly.
 *
 * @see FRTShellS3
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell_s3.h"

namespace fem::model {

/**
 * Constructs the three-node MITC3 finite-rotation shell.
 *
 * The element uses linear triangular geometry and cubic area quadrature. The
 * oriented surface representations are created on demand by `surface()`.
 *
 * @param id Element identifier.
 * @param nodes Three nodal identifiers in positive triangular ordering.
 */
FRTShellS3::FRTShellS3(ID id, const std::array<ID, 3>& nodes)
    : FRTShell<3>       (id, nodes),
      integration_scheme_(math::quadrature::Domain::DOMAIN_ISO_TRI,
                          math::quadrature::Order::ORDER_CUBIC) {}

/**
 * Returns the public element type identifier.
 *
 * @return String identifier `MITC3FRT`.
 */
std::string FRTShellS3::type_name() const {
    return "MITC3FRT";
}

/**
 * Returns an oriented triangular surface representation.
 *
 * Surface one follows the positive shell orientation. Any opposite surface
 * identifier reverses the triangle ordering while preserving the first node.
 *
 * @param surface_id Requested shell side identifier.
 * @return Oriented three-node surface object.
 */
std::shared_ptr<SurfaceInterface> FRTShellS3::surface(int surface_id) {
    return std::make_shared<Surface3>(
        surface_id == 1
            ? std::array<ID, 3>{this->nodes()[0], this->nodes()[1], this->nodes()[2]}
            : std::array<ID, 3>{this->nodes()[0], this->nodes()[2], this->nodes()[1]}
    );
}

/**
 * Returns the numerical area quadrature used by the element.
 *
 * @return Cubic triangular quadrature rule.
 */
const math::quadrature::Quadrature& FRTShellS3::integration_scheme() const {
    return integration_scheme_;
}

/**
 * Returns all numerical integration-point coordinates for result storage.
 *
 * The first two columns contain the natural triangle coordinates. The third
 * through-thickness coordinate is zero because generalized shell resultants
 * belong to the midsurface.
 *
 * @return Matrix of `(r,s,t)` integration-point coordinates.
 */
RowMatrix FRTShellS3::stress_strain_ip_rst() {
    RowMatrix rst(integration_scheme_.count(), 3);
    rst.setZero();

    // Copy the natural coordinates from the quadrature rule in its native order
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);
        rst(ip, 0) = point.r;
        rst(ip, 1) = point.s;
    }

    return rst;
}

/**
 * Returns the three natural nodal coordinates for nodal result recovery.
 *
 * @return Matrix containing the three triangle vertices in `(r,s,t)` format.
 */
RowMatrix FRTShellS3::stress_strain_nodal_rst() {
    RowMatrix rst(3, 3);
    rst << Precision(0), Precision(0), Precision(0),
           Precision(1), Precision(0), Precision(0),
           Precision(0), Precision(1), Precision(0);
    return rst;
}

/**
 * Evaluates the linear barycentric shape functions on the reference triangle.
 *
 * @param r First natural triangle coordinate.
 * @param s Second natural triangle coordinate.
 * @return Three linear shape-function values.
 */
FRTShellS3::VecN FRTShellS3::shape_function(Precision r, Precision s) const {
    VecN shape;
    shape << Precision(1) - r - s,
             r,
             s;
    return shape;
}

/**
 * Evaluates the constant natural derivatives of the linear shape functions.
 *
 * @param r Unused first natural coordinate because the derivatives are constant.
 * @param s Unused second natural coordinate because the derivatives are constant.
 * @return Matrix whose rows contain `[dN_i/dr,dN_i/ds]`.
 */
FRTShellS3::MatN2 FRTShellS3::shape_derivative(Precision r, Precision s) const {
    (void) r;
    (void) s;

    MatN2 derivative;
    derivative << Precision(-1), Precision(-1),
                  Precision( 1), Precision( 0),
                  Precision( 0), Precision( 1);
    return derivative;
}

/**
 * Returns the natural coordinates of the three triangle vertices.
 *
 * @return Matrix containing `(0,0)`, `(1,0)` and `(0,1)`.
 */
FRTShellS3::MatN2 FRTShellS3::node_coords_natural() const {
    MatN2 coordinates;
    coordinates << Precision(0), Precision(0),
                   Precision(1), Precision(0),
                   Precision(0), Precision(1);
    return coordinates;
}

/**
 * Returns the three edge-midpoint tying positions of the MITC3 shear field.
 *
 * The ordering is
 *
 * - 0: bottom edge `(1/2,0)`, tying `gamma_r3`,
 * - 1: left edge `(0,1/2)`, tying `gamma_s3`,
 * - 2: hypotenuse `(1/2,1/2)`, tying the third edge direction.
 *
 * @return Ordered natural tying-point coordinates.
 */
std::vector<Vec2> FRTShellS3::tying_point_coordinates() const {
    return {
        Vec2(Precision(0.5), Precision(0.0)),
        Vec2(Precision(0.0), Precision(0.5)),
        Vec2(Precision(0.5), Precision(0.5))
    };
}

/**
 * Applies the isotropic MITC3 transverse-shear interpolation.
 *
 * With the three edge-midpoint tying values, the assumed natural shear field is
 *
 *     gamma_r3 = gamma_r3^(1) + c s,
 *     gamma_s3 = gamma_s3^(2) - c r,
 *
 * where
 *
 *     c = gamma_s3^(2) - gamma_r3^(1)
 *       - gamma_s3^(3) + gamma_r3^(3).
 *
 * The same coefficient is applied to the B rows, preserving exact consistency
 * between the assumed strain and its first derivative.
 *
 * @param data Active tying-point strain and B data.
 * @param point Target natural shell point.
 * @param strain_nat Compatible natural strain vector to modify.
 * @param B_nat Optional compatible natural B matrix to modify.
 */
void FRTShellS3::apply_mitc_natural(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat
) const {
    constexpr Index gamma_r3 = 6;
    constexpr Index gamma_s3 = 7;

    const Precision r = point.r;
    const Precision s = point.s;

    // Build the coefficient of the common isotropic linear shear variation
    const Precision c = data.tying_strain_nat[1](gamma_s3)
                      - data.tying_strain_nat[0](gamma_r3)
                      - data.tying_strain_nat[2](gamma_s3)
                      + data.tying_strain_nat[2](gamma_r3);

    strain_nat(gamma_r3) = data.tying_strain_nat[0](gamma_r3) + s * c;
    strain_nat(gamma_s3) = data.tying_strain_nat[1](gamma_s3) - r * c;

    if (B_nat) {
        // Apply the identical linear interpolation to the two shear B rows
        const StaticMatrix<1, num_dofs> c_B =
              data.tying_B_nat[1].row(gamma_s3)
            - data.tying_B_nat[0].row(gamma_r3)
            - data.tying_B_nat[2].row(gamma_s3)
            + data.tying_B_nat[2].row(gamma_r3);

        B_nat->row(gamma_r3) = data.tying_B_nat[0].row(gamma_r3) + s * c_B;
        B_nat->row(gamma_s3) = data.tying_B_nat[1].row(gamma_s3) - r * c_B;
    }
}

/**
 * Applies the transpose of the MITC3 assumed-shear interpolation.
 *
 * Membrane and curvature weights remain compatible at the integration point.
 * The two assumed shear weights are distributed to the three edge-midpoint
 * compatible shear Hessians using the exact transpose of
 * `apply_mitc_natural()`.
 *
 * @param point Integration point defining the interpolation coefficients.
 * @param assumed_weights Generalized natural weights after local-basis pull-back.
 * @param compatible_weights Compatible integration-point weights to increment.
 * @param tying_weights Zeroed compatible tying-point weights to increment.
 */
void FRTShellS3::pull_back_mitc_resultants(
    const ReferencePoint& point,
    const Vec8&           assumed_weights,
    Vec8&                 compatible_weights,
    Span<Vec8>            tying_weights
) const {
    constexpr Index gamma_r3 = 6;
    constexpr Index gamma_s3 = 7;

    const Precision r  = point.r;
    const Precision s  = point.s;
    const Precision wr = assumed_weights(gamma_r3);
    const Precision ws = assumed_weights(gamma_s3);

    // The first six compatible components are not replaced by MITC3
    compatible_weights.template segment<6>(0) +=
        assumed_weights.template segment<6>(0);

    // Distribute the assumed shear weights through the transposed interpolation
    tying_weights[0](gamma_r3) += (Precision(1) - s) * wr + r * ws;
    tying_weights[1](gamma_s3) += s * wr + (Precision(1) - r) * ws;
    tying_weights[2](gamma_r3) += s * wr - r * ws;
    tying_weights[2](gamma_s3) += -s * wr + r * ws;
}

} // namespace fem::model
