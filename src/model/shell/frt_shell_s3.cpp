/**
 * @file frt_shell_s3.cpp
 * @brief Implements the three-node finite-rotation MITC3 shell element.
 *
 * The assumed covariant shear field follows the isotropic MITC3 interpolation
 * on the natural triangle. The same coefficients are applied to the shear
 * values and their B rows; geometric-tangent weights are pulled back through
 * the transposed interpolation.
 *
 * @see FRTShellS3
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell_s3.h"

namespace fem::model {

FRTShellS3::FRTShellS3(ID id, const std::array<ID, 3>& nodes)
    : FRTShell<3>       (id, nodes),
      geometry          (nodes),
      integration_scheme_(math::quadrature::Domain::DOMAIN_ISO_TRI,
                          math::quadrature::Order::ORDER_CUBIC) {}

std::string FRTShellS3::type_name() const {
    return "MITC3FRT";
}

std::shared_ptr<SurfaceInterface> FRTShellS3::surface(int surface_id) {
    return std::make_shared<Surface3>(
        surface_id == 1
            ? std::array<ID, 3>{this->nodes()[0], this->nodes()[1], this->nodes()[2]}
            : std::array<ID, 3>{this->nodes()[0], this->nodes()[2], this->nodes()[1]}
    );
}

const math::quadrature::Quadrature& FRTShellS3::integration_scheme() const {
    return integration_scheme_;
}

RowMatrix FRTShellS3::stress_strain_ip_rst() {
    RowMatrix rst(integration_scheme_.count(), 3);
    rst.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);
        rst(ip, 0) = point.r;
        rst(ip, 1) = point.s;
    }

    return rst;
}

RowMatrix FRTShellS3::stress_strain_nodal_rst() {
    RowMatrix rst(3, 3);
    rst << Precision(0), Precision(0), Precision(0),
           Precision(1), Precision(0), Precision(0),
           Precision(0), Precision(1), Precision(0);
    return rst;
}

/**
 * Evaluates the linear barycentric shape functions on the reference triangle.
 */
FRTShellS3::VecN FRTShellS3::shape_function(Precision r, Precision s) const {
    VecN shape;
    shape << Precision(1) - r - s,
             r,
             s;
    return shape;
}

/**
 * Returns the constant natural derivatives of the linear triangle.
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
 * The ordering is:
 *
 * - 0: bottom edge, `(1/2,0)`, tying `gamma_r3`,
 * - 1: left edge, `(0,1/2)`, tying `gamma_s3`,
 * - 2: hypotenuse, `(1/2,1/2)`, tying the third edge direction.
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
 * This form gives identical linear variation in all three edge directions.
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
        const StaticMatrix<1, num_dofs> c_B =
              data.tying_B_nat[1].row(gamma_s3)
            - data.tying_B_nat[0].row(gamma_r3)
            - data.tying_B_nat[2].row(gamma_s3)
            + data.tying_B_nat[2].row(gamma_r3);

        B_nat->row(gamma_r3) = data.tying_B_nat[0].row(gamma_r3) + s * c_B;
        B_nat->row(gamma_s3) = data.tying_B_nat[1].row(gamma_s3) - r * c_B;
    }
}

void FRTShellS3::pull_back_mitc_resultants(
    const ReferencePoint& point,
    const Vec8&           assumed_weights,
    Vec8&                 compatible_weights,
    std::vector<Vec8>&    tying_weights
) const {
    constexpr Index gamma_r3 = 6;
    constexpr Index gamma_s3 = 7;

    const Precision r = point.r;
    const Precision s = point.s;
    const Precision wr = assumed_weights(gamma_r3);
    const Precision ws = assumed_weights(gamma_s3);

    compatible_weights.template segment<6>(0) += assumed_weights.template segment<6>(0);

    tying_weights[0](gamma_r3) += (Precision(1) - s) * wr + r * ws;
    tying_weights[1](gamma_s3) += s * wr + (Precision(1) - r) * ws;
    tying_weights[2](gamma_r3) += s * wr - r * ws;
    tying_weights[2](gamma_s3) += -s * wr + r * ws;
}

} // namespace fem::model
