/**
 * @file frt_shell_s4.cpp
 * @brief Implements the four-node finite-rotation MITC4 shell element.
 *
 * The MITC4 interpolation is performed in natural covariant components before
 * the common shell kernel transforms the strains into the pointwise local
 * reference basis. Values, B rows and G matrices use exactly the same tying
 * coefficients.
 *
 * @see FRTShellS4
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell_s4.h"

namespace fem::model {

FRTShellS4::FRTShellS4(ID id, const std::array<ID, 4>& nodes)
    : FRTShell<4>       (id, nodes),
      geometry          (nodes),
      integration_scheme_(math::quadrature::Domain::DOMAIN_ISO_QUAD,
                          math::quadrature::Order::ORDER_CUBIC) {}

std::string FRTShellS4::type_name() const {
    return "MITC4FRT";
}

std::shared_ptr<SurfaceInterface> FRTShellS4::surface(int surface_id) {
    return std::make_shared<Surface4>(
        surface_id == 1
            ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
            : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]}
    );
}

const math::quadrature::Quadrature& FRTShellS4::integration_scheme() const {
    return integration_scheme_;
}

RowMatrix FRTShellS4::stress_strain_ip_rst() {
    RowMatrix rst(integration_scheme_.count(), 3);
    rst.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);
        rst(ip, 0) = point.r;
        rst(ip, 1) = point.s;
    }

    return rst;
}

RowMatrix FRTShellS4::stress_strain_nodal_rst() {
    RowMatrix rst(4, 3);
    rst << Precision(-1), Precision(-1), Precision(0),
           Precision( 1), Precision(-1), Precision(0),
           Precision( 1), Precision( 1), Precision(0),
           Precision(-1), Precision( 1), Precision(0);
    return rst;
}

FRTShellS4::VecN FRTShellS4::shape_function(Precision r, Precision s) const {
    VecN shape;
    shape(0) = Precision(0.25) * (Precision(1) - r) * (Precision(1) - s);
    shape(1) = Precision(0.25) * (Precision(1) + r) * (Precision(1) - s);
    shape(2) = Precision(0.25) * (Precision(1) + r) * (Precision(1) + s);
    shape(3) = Precision(0.25) * (Precision(1) - r) * (Precision(1) + s);
    return shape;
}

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
 * `gamma_r3` is interpolated between the bottom and top tying points using the
 * natural coordinate `s`. `gamma_s3` is interpolated between the left and right
 * tying points using `r`.
 */
void FRTShellS4::apply_mitc_natural(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat,
    Vec6NMat*             G_nat
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
        B_nat->row(gamma_r3) = bottom * data.tying_B_nat[0].row(gamma_r3)
                             + top    * data.tying_B_nat[1].row(gamma_r3);
        B_nat->row(gamma_s3) = left   * data.tying_B_nat[2].row(gamma_s3)
                             + right  * data.tying_B_nat[3].row(gamma_s3);
    }

    if (G_nat) {
        (*G_nat)[gamma_r3] = bottom * data.tying_G_nat[0][gamma_r3]
                           + top    * data.tying_G_nat[1][gamma_r3];
        (*G_nat)[gamma_s3] = left   * data.tying_G_nat[2][gamma_s3]
                           + right  * data.tying_G_nat[3][gamma_s3];
    }
}

} // namespace fem::model
