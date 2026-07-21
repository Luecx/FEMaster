/**
 * @file frt_shell_s6.cpp
 * @brief Implements the six-node finite-rotation MITC6-b shell element.
 *
 * The assumed in-plane strain field follows the isotropic triangular
 * interpolation of Lee and Bathe. The transverse shear field uses the MITC6-b
 * linear interpolation with two tying points on each edge. Membrane strain and
 * curvature are interpolated independently with the same in-plane operator.
 *
 * @see FRTShellS6
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell_s6.h"

#include <cmath>

namespace fem::model {

FRTShellS6::FRTShellS6(ID id, const std::array<ID, 6>& nodes)
    : FRTShell<6>       (id, nodes),
      geometry          (nodes),
      integration_scheme_(math::quadrature::Domain::DOMAIN_ISO_TRI,
                          math::quadrature::Order::ORDER_CUBIC) {}

std::string FRTShellS6::type_name() const {
    return "MITC6FRT";
}

std::shared_ptr<SurfaceInterface> FRTShellS6::surface(int surface_id) {
    return std::make_shared<Surface6>(
        surface_id == 1
            ? std::array<ID, 6>{this->nodes()[0], this->nodes()[1], this->nodes()[2],
                                this->nodes()[3], this->nodes()[4], this->nodes()[5]}
            : std::array<ID, 6>{this->nodes()[0], this->nodes()[2], this->nodes()[1],
                                this->nodes()[5], this->nodes()[4], this->nodes()[3]}
    );
}

const math::quadrature::Quadrature& FRTShellS6::integration_scheme() const {
    return integration_scheme_;
}

RowMatrix FRTShellS6::stress_strain_ip_rst() {
    RowMatrix rst(integration_scheme_.count(), 3);
    rst.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);
        rst(ip, 0) = point.r;
        rst(ip, 1) = point.s;
    }

    return rst;
}

RowMatrix FRTShellS6::stress_strain_nodal_rst() {
    RowMatrix rst(6, 3);
    rst << Precision(0.0), Precision(0.0), Precision(0),
           Precision(1.0), Precision(0.0), Precision(0),
           Precision(0.0), Precision(1.0), Precision(0),
           Precision(0.5), Precision(0.0), Precision(0),
           Precision(0.5), Precision(0.5), Precision(0),
           Precision(0.0), Precision(0.5), Precision(0);
    return rst;
}

FRTShellS6::VecN FRTShellS6::shape_function(Precision r, Precision s) const {
    const Precision t = Precision(1) - r - s;

    VecN shape;
    shape << t * (Precision(2) * t - Precision(1)),
             r * (Precision(2) * r - Precision(1)),
             s * (Precision(2) * s - Precision(1)),
             Precision(4) * r * t,
             Precision(4) * r * s,
             Precision(4) * s * t;
    return shape;
}

FRTShellS6::MatN2 FRTShellS6::shape_derivative(Precision r, Precision s) const {
    MatN2 derivative;

    derivative(0, 0) = -Precision(3) + Precision(4) * (r + s);
    derivative(0, 1) = -Precision(3) + Precision(4) * (r + s);
    derivative(1, 0) =  Precision(4) * r - Precision(1);
    derivative(1, 1) =  Precision(0);
    derivative(2, 0) =  Precision(0);
    derivative(2, 1) =  Precision(4) * s - Precision(1);
    derivative(3, 0) =  Precision(4) - Precision(8) * r - Precision(4) * s;
    derivative(3, 1) = -Precision(4) * r;
    derivative(4, 0) =  Precision(4) * s;
    derivative(4, 1) =  Precision(4) * r;
    derivative(5, 0) = -Precision(4) * s;
    derivative(5, 1) =  Precision(4) - Precision(4) * r - Precision(8) * s;

    return derivative;
}

FRTShellS6::MatN2 FRTShellS6::node_coords_natural() const {
    MatN2 coordinates;
    coordinates << Precision(0.0), Precision(0.0),
                   Precision(1.0), Precision(0.0),
                   Precision(0.0), Precision(1.0),
                   Precision(0.5), Precision(0.0),
                   Precision(0.5), Precision(0.5),
                   Precision(0.0), Precision(0.5);
    return coordinates;
}

/**
 * Returns the MITC6-b in-plane and transverse-shear tying positions.
 *
 * The constants
 *
 *     r1 = s1 = 1/2 - 1/(2 sqrt(3)),
 *     r2 = s2 = 1/2 + 1/(2 sqrt(3))
 *
 * are the two-point Gauss positions mapped onto every triangle edge. The first
 * nine points define the isotropic in-plane strain interpolation. The final six
 * points define the MITC6-b transverse-shear field.
 */
std::vector<Vec2> FRTShellS6::tying_point_coordinates() const {
    const Precision sqrt3     = std::sqrt(Precision(3));
    const Precision r1        = Precision(0.5) - Precision(0.5) / sqrt3;
    const Precision r2        = Precision(0.5) + Precision(0.5) / sqrt3;
    const Precision inv_sqrt3 = Precision(1) / sqrt3;

    return {
        // In-plane e_rr sampling
        Vec2(r1, Precision(0)),
        Vec2(r2, Precision(0)),
        Vec2(r1, inv_sqrt3),

        // In-plane e_ss sampling
        Vec2(Precision(0), r1),
        Vec2(Precision(0), r2),
        Vec2(inv_sqrt3, r1),

        // In-plane e_qq sampling on and inside the hypotenuse
        Vec2(r2, r1),
        Vec2(r1, r2),
        Vec2(r1, r1),

        // Transverse shear on the r, s and q edges
        Vec2(r1, Precision(0)),
        Vec2(r2, Precision(0)),
        Vec2(Precision(0), r1),
        Vec2(Precision(0), r2),
        Vec2(r2, r1),
        Vec2(r1, r2)
    };
}

/**
 * Applies the MITC6-b in-plane and transverse-shear interpolation.
 *
 * For both membrane strain and curvature, three normal strain fields
 * `e_rr`, `e_ss` and `e_qq` are interpolated isotropically. The engineering
 * shear follows from `gamma_rs = e_rr + e_ss - 2 e_qq`.
 *
 * The transverse-shear field is linear inside the triangle and uses two tying
 * values on each edge. Values, B rows and G matrices follow exactly the same
 * coefficient construction.
 */
void FRTShellS6::apply_mitc_natural(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat,
    Vec6NMat*             G_nat
) const {
    const Precision sqrt3 = std::sqrt(Precision(3));
    const Precision r1    = Precision(0.5) - Precision(0.5) / sqrt3;
    const Precision r     = point.r;
    const Precision s     = point.s;

    // Apply the isotropic in-plane interpolation independently to membrane
    // strain and curvature
    auto interpolate_in_plane_values = [&](Index start) {
        const auto q_value = [&](Index tying) {
            const Vec8& value = data.tying_strain_nat[static_cast<std::size_t>(tying)];
            return Precision(0.5)
                 * (value(start + 0) + value(start + 1) - value(start + 2));
        };

        const Precision m_rr = Precision(0.5)
                             * (data.tying_strain_nat[0](start + 0)
                              + data.tying_strain_nat[1](start + 0));
        const Precision l_rr = Precision(0.5) * sqrt3
                             * (data.tying_strain_nat[1](start + 0)
                              - data.tying_strain_nat[0](start + 0));

        const Precision m_ss = Precision(0.5)
                             * (data.tying_strain_nat[3](start + 1)
                              + data.tying_strain_nat[4](start + 1));
        const Precision l_ss = Precision(0.5) * sqrt3
                             * (data.tying_strain_nat[4](start + 1)
                              - data.tying_strain_nat[3](start + 1));

        const Precision m_qq = Precision(0.5) * (q_value(6) + q_value(7));
        const Precision l_qq = Precision(0.5) * sqrt3 * (q_value(7) - q_value(6));

        const Precision a1 = m_rr - l_rr;
        const Precision b1 = Precision(2) * l_rr;
        const Precision c1 = sqrt3
                           * (data.tying_strain_nat[2](start + 0) - a1 - b1 * r1);

        const Precision a2 = m_ss - l_ss;
        const Precision c2 = Precision(2) * l_ss;
        const Precision b2 = sqrt3
                           * (data.tying_strain_nat[5](start + 1) - a2 - c2 * r1);

        const Precision a3 = m_qq + l_qq;
        const Precision b3 = -Precision(2) * l_qq;
        const Precision c3 = sqrt3 * (q_value(8) - a3 - b3 * r1);

        const Precision e_rr = a1 + b1 * r + c1 * s;
        const Precision e_ss = a2 + b2 * r + c2 * s;
        const Precision e_qq = a3 + b3 * r + c3 * (Precision(1) - r - s);

        strain_nat(start + 0) = e_rr;
        strain_nat(start + 1) = e_ss;
        strain_nat(start + 2) = e_rr + e_ss - Precision(2) * e_qq;
    };

    interpolate_in_plane_values(0);
    interpolate_in_plane_values(3);

    if (B_nat) {
        auto interpolate_in_plane_B = [&](Index start) {
            const auto q_row = [&](Index tying) -> StaticMatrix<1, num_dofs> {
                StaticMatrix<1, num_dofs> result =
                      Precision(0.5) * data.tying_B_nat[tying].row(start + 0)
                    + Precision(0.5) * data.tying_B_nat[tying].row(start + 1)
                    - Precision(0.5) * data.tying_B_nat[tying].row(start + 2);
                return result;
            };

            const StaticMatrix<1, num_dofs> m_rr =
                  Precision(0.5) * data.tying_B_nat[0].row(start + 0)
                + Precision(0.5) * data.tying_B_nat[1].row(start + 0);
            const StaticMatrix<1, num_dofs> l_rr =
                  Precision(0.5) * sqrt3 * data.tying_B_nat[1].row(start + 0)
                - Precision(0.5) * sqrt3 * data.tying_B_nat[0].row(start + 0);

            const StaticMatrix<1, num_dofs> m_ss =
                  Precision(0.5) * data.tying_B_nat[3].row(start + 1)
                + Precision(0.5) * data.tying_B_nat[4].row(start + 1);
            const StaticMatrix<1, num_dofs> l_ss =
                  Precision(0.5) * sqrt3 * data.tying_B_nat[4].row(start + 1)
                - Precision(0.5) * sqrt3 * data.tying_B_nat[3].row(start + 1);

            const StaticMatrix<1, num_dofs> q_6 = q_row(6);
            const StaticMatrix<1, num_dofs> q_7 = q_row(7);
            const StaticMatrix<1, num_dofs> q_8 = q_row(8);

            const StaticMatrix<1, num_dofs> m_qq = Precision(0.5) * (q_6 + q_7);
            const StaticMatrix<1, num_dofs> l_qq = Precision(0.5) * sqrt3 * (q_7 - q_6);

            const StaticMatrix<1, num_dofs> a1 = m_rr - l_rr;
            const StaticMatrix<1, num_dofs> b1 = Precision(2) * l_rr;
            const StaticMatrix<1, num_dofs> c1 =
                sqrt3 * (data.tying_B_nat[2].row(start + 0) - a1 - b1 * r1);

            const StaticMatrix<1, num_dofs> a2 = m_ss - l_ss;
            const StaticMatrix<1, num_dofs> c2 = Precision(2) * l_ss;
            const StaticMatrix<1, num_dofs> b2 =
                sqrt3 * (data.tying_B_nat[5].row(start + 1) - a2 - c2 * r1);

            const StaticMatrix<1, num_dofs> a3 = m_qq + l_qq;
            const StaticMatrix<1, num_dofs> b3 = -Precision(2) * l_qq;
            const StaticMatrix<1, num_dofs> c3 = sqrt3 * (q_8 - a3 - b3 * r1);

            const StaticMatrix<1, num_dofs> e_rr = a1 + b1 * r + c1 * s;
            const StaticMatrix<1, num_dofs> e_ss = a2 + b2 * r + c2 * s;
            const StaticMatrix<1, num_dofs> e_qq =
                a3 + b3 * r + c3 * (Precision(1) - r - s);

            B_nat->row(start + 0) = e_rr;
            B_nat->row(start + 1) = e_ss;
            B_nat->row(start + 2) = e_rr + e_ss - Precision(2) * e_qq;
        };

        interpolate_in_plane_B(0);
        interpolate_in_plane_B(3);
    }

    if (G_nat) {
        auto interpolate_in_plane_G = [&](Index start) {
            const auto q_matrix = [&](Index tying) {
                const Mat6N result = Precision(0.5)
                    * (data.tying_G_nat[tying][start + 0]
                     + data.tying_G_nat[tying][start + 1]
                     - data.tying_G_nat[tying][start + 2]);
                return result;
            };

            const Mat6N m_rr = Precision(0.5)
                             * (data.tying_G_nat[0][start + 0]
                              + data.tying_G_nat[1][start + 0]);
            const Mat6N l_rr = Precision(0.5) * sqrt3
                             * (data.tying_G_nat[1][start + 0]
                              - data.tying_G_nat[0][start + 0]);

            const Mat6N m_ss = Precision(0.5)
                             * (data.tying_G_nat[3][start + 1]
                              + data.tying_G_nat[4][start + 1]);
            const Mat6N l_ss = Precision(0.5) * sqrt3
                             * (data.tying_G_nat[4][start + 1]
                              - data.tying_G_nat[3][start + 1]);

            const Mat6N m_qq = Precision(0.5) * (q_matrix(6) + q_matrix(7));
            const Mat6N l_qq = Precision(0.5) * sqrt3 * (q_matrix(7) - q_matrix(6));

            const Mat6N a1 = m_rr - l_rr;
            const Mat6N b1 = Precision(2) * l_rr;
            const Mat6N c1 = sqrt3
                           * (data.tying_G_nat[2][start + 0] - a1 - b1 * r1);

            const Mat6N a2 = m_ss - l_ss;
            const Mat6N c2 = Precision(2) * l_ss;
            const Mat6N b2 = sqrt3
                           * (data.tying_G_nat[5][start + 1] - a2 - c2 * r1);

            const Mat6N a3 = m_qq + l_qq;
            const Mat6N b3 = -Precision(2) * l_qq;
            const Mat6N c3 = sqrt3 * (q_matrix(8) - a3 - b3 * r1);

            const Mat6N e_rr = a1 + b1 * r + c1 * s;
            const Mat6N e_ss = a2 + b2 * r + c2 * s;
            const Mat6N e_qq = a3 + b3 * r + c3 * (Precision(1) - r - s);

            (*G_nat)[start + 0] = e_rr;
            (*G_nat)[start + 1] = e_ss;
            (*G_nat)[start + 2] = e_rr + e_ss - Precision(2) * e_qq;
        };

        interpolate_in_plane_G(0);
        interpolate_in_plane_G(3);
    }

    // MITC6-b transverse-shear interpolation
    constexpr Index gamma_r3 = 6;
    constexpr Index gamma_s3 = 7;

    const auto mean = [](Precision a, Precision b) {
        return Precision(0.5) * (a + b);
    };
    const auto slope = [sqrt3](Precision a, Precision b) {
        return Precision(0.5) * sqrt3 * (b - a);
    };

    const Precision m_rt_1 = mean(data.tying_strain_nat[9](gamma_r3),
                                  data.tying_strain_nat[10](gamma_r3));
    const Precision l_rt_1 = slope(data.tying_strain_nat[9](gamma_r3),
                                   data.tying_strain_nat[10](gamma_r3));
    const Precision m_st_2 = mean(data.tying_strain_nat[11](gamma_s3),
                                  data.tying_strain_nat[12](gamma_s3));
    const Precision l_st_2 = slope(data.tying_strain_nat[11](gamma_s3),
                                   data.tying_strain_nat[12](gamma_s3));

    const Precision m_rt_3 = mean(data.tying_strain_nat[13](gamma_r3),
                                  data.tying_strain_nat[14](gamma_r3));
    const Precision l_rt_3 = slope(data.tying_strain_nat[13](gamma_r3),
                                   data.tying_strain_nat[14](gamma_r3));
    const Precision m_st_3 = mean(data.tying_strain_nat[13](gamma_s3),
                                  data.tying_strain_nat[14](gamma_s3));
    const Precision l_st_3 = slope(data.tying_strain_nat[13](gamma_s3),
                                   data.tying_strain_nat[14](gamma_s3));

    const Precision a1 = m_rt_1 - l_rt_1;
    const Precision b1 = Precision(2) * l_rt_1;
    const Precision a2 = m_st_2 - l_st_2;
    const Precision c2 = Precision(2) * l_st_2;
    const Precision c1 = (a2 + c2 - a1) - (m_st_3 + l_st_3 - m_rt_3 - l_rt_3);
    const Precision b2 = (a1 + b1 - a2) + (m_st_3 - l_st_3 - m_rt_3 + l_rt_3);

    strain_nat(gamma_r3) = a1 + b1 * r + c1 * s;
    strain_nat(gamma_s3) = a2 + b2 * r + c2 * s;

    if (B_nat) {
        const auto mean_row = [](const StaticMatrix<1, num_dofs>& a,
                                 const StaticMatrix<1, num_dofs>& b) {
            return StaticMatrix<1, num_dofs>(Precision(0.5) * (a + b));
        };
        const auto slope_row = [sqrt3](const StaticMatrix<1, num_dofs>& a,
                                       const StaticMatrix<1, num_dofs>& b) {
            return StaticMatrix<1, num_dofs>(Precision(0.5) * sqrt3 * (b - a));
        };

        const StaticMatrix<1, num_dofs> m_rt1 =
            mean_row(data.tying_B_nat[9].row(gamma_r3),
                     data.tying_B_nat[10].row(gamma_r3));
        const StaticMatrix<1, num_dofs> l_rt1 =
            slope_row(data.tying_B_nat[9].row(gamma_r3),
                      data.tying_B_nat[10].row(gamma_r3));
        const StaticMatrix<1, num_dofs> m_st2 =
            mean_row(data.tying_B_nat[11].row(gamma_s3),
                     data.tying_B_nat[12].row(gamma_s3));
        const StaticMatrix<1, num_dofs> l_st2 =
            slope_row(data.tying_B_nat[11].row(gamma_s3),
                      data.tying_B_nat[12].row(gamma_s3));
        const StaticMatrix<1, num_dofs> m_rt3 =
            mean_row(data.tying_B_nat[13].row(gamma_r3),
                     data.tying_B_nat[14].row(gamma_r3));
        const StaticMatrix<1, num_dofs> l_rt3 =
            slope_row(data.tying_B_nat[13].row(gamma_r3),
                      data.tying_B_nat[14].row(gamma_r3));
        const StaticMatrix<1, num_dofs> m_st3 =
            mean_row(data.tying_B_nat[13].row(gamma_s3),
                     data.tying_B_nat[14].row(gamma_s3));
        const StaticMatrix<1, num_dofs> l_st3 =
            slope_row(data.tying_B_nat[13].row(gamma_s3),
                      data.tying_B_nat[14].row(gamma_s3));

        const StaticMatrix<1, num_dofs> A1 = m_rt1 - l_rt1;
        const StaticMatrix<1, num_dofs> B1 = Precision(2) * l_rt1;
        const StaticMatrix<1, num_dofs> A2 = m_st2 - l_st2;
        const StaticMatrix<1, num_dofs> C2 = Precision(2) * l_st2;
        const StaticMatrix<1, num_dofs> C1 =
            (A2 + C2 - A1) - (m_st3 + l_st3 - m_rt3 - l_rt3);
        const StaticMatrix<1, num_dofs> B2 =
            (A1 + B1 - A2) + (m_st3 - l_st3 - m_rt3 + l_rt3);

        B_nat->row(gamma_r3) = A1 + B1 * r + C1 * s;
        B_nat->row(gamma_s3) = A2 + B2 * r + C2 * s;
    }

    if (G_nat) {
        const auto mean_matrix = [](const Mat6N& a, const Mat6N& b) {
            const Mat6N result = Precision(0.5) * (a + b);
            return result;
        };
        const auto slope_matrix = [sqrt3](const Mat6N& a, const Mat6N& b) {
            const Mat6N result = Precision(0.5) * sqrt3 * (b - a);
            return result;
        };

        const Mat6N m_rt1 = mean_matrix(data.tying_G_nat[9][gamma_r3],
                                        data.tying_G_nat[10][gamma_r3]);
        const Mat6N l_rt1 = slope_matrix(data.tying_G_nat[9][gamma_r3],
                                         data.tying_G_nat[10][gamma_r3]);
        const Mat6N m_st2 = mean_matrix(data.tying_G_nat[11][gamma_s3],
                                        data.tying_G_nat[12][gamma_s3]);
        const Mat6N l_st2 = slope_matrix(data.tying_G_nat[11][gamma_s3],
                                         data.tying_G_nat[12][gamma_s3]);
        const Mat6N m_rt3 = mean_matrix(data.tying_G_nat[13][gamma_r3],
                                        data.tying_G_nat[14][gamma_r3]);
        const Mat6N l_rt3 = slope_matrix(data.tying_G_nat[13][gamma_r3],
                                         data.tying_G_nat[14][gamma_r3]);
        const Mat6N m_st3 = mean_matrix(data.tying_G_nat[13][gamma_s3],
                                        data.tying_G_nat[14][gamma_s3]);
        const Mat6N l_st3 = slope_matrix(data.tying_G_nat[13][gamma_s3],
                                         data.tying_G_nat[14][gamma_s3]);

        const Mat6N A1 = m_rt1 - l_rt1;
        const Mat6N B1 = Precision(2) * l_rt1;
        const Mat6N A2 = m_st2 - l_st2;
        const Mat6N C2 = Precision(2) * l_st2;
        const Mat6N C1 = (A2 + C2 - A1) - (m_st3 + l_st3 - m_rt3 - l_rt3);
        const Mat6N B2 = (A1 + B1 - A2) + (m_st3 - l_st3 - m_rt3 + l_rt3);

        (*G_nat)[gamma_r3] = A1 + B1 * r + C1 * s;
        (*G_nat)[gamma_s3] = A2 + B2 * r + C2 * s;
    }
}

} // namespace fem::model
