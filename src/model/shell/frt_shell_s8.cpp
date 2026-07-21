/**
 * @file frt_shell_s8.cpp
 * @brief Implements the eight-node finite-rotation MITC8 shell element.
 *
 * The implementation follows the classical MITC8 construction. Eight
 * in-layer sampling tensors are interpolated with serendipity functions on the
 * reduced square `[-a,a] x [-a,a]`, where `a = 1/sqrt(3)`. At the midside
 * sampling points, one directly evaluated normal component is retained while
 * the remaining tensor components are obtained from the adjacent corner
 * tensors in an orthogonalized reference-surface basis.
 *
 * The two covariant transverse-shear fields use four boundary sampling values
 * and one internal function. The fifth coefficient is the mean of the two
 * associated interior sampling values. Every assumed-strain interpolation is
 * applied identically to the strain values, B rows and G matrices.
 *
 * @see FRTShellS8
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell_s8.h"

#include "../../core/logging.h"

#include <cmath>

namespace fem::model {

FRTShellS8::FRTShellS8(ID id, const std::array<ID, 8>& nodes)
    : FRTShell<8>       (id, nodes),
      geometry          (nodes),
      integration_scheme_(math::quadrature::Domain::DOMAIN_ISO_QUAD,
                          math::quadrature::Order::ORDER_QUINTIC) {}

std::string FRTShellS8::type_name() const {
    return "MITC8FRT";
}

std::shared_ptr<SurfaceInterface> FRTShellS8::surface(int surface_id) {
    return std::make_shared<Surface8>(
        surface_id == 1
            ? std::array<ID, 8>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3],
                                this->nodes()[4], this->nodes()[5], this->nodes()[6], this->nodes()[7]}
            : std::array<ID, 8>{this->nodes()[0], this->nodes()[3], this->nodes()[2], this->nodes()[1],
                                this->nodes()[7], this->nodes()[6], this->nodes()[5], this->nodes()[4]}
    );
}

const math::quadrature::Quadrature& FRTShellS8::integration_scheme() const {
    return integration_scheme_;
}

RowMatrix FRTShellS8::stress_strain_ip_rst() {
    RowMatrix rst(integration_scheme_.count(), 3);
    rst.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);
        rst(ip, 0) = point.r;
        rst(ip, 1) = point.s;
    }

    return rst;
}

RowMatrix FRTShellS8::stress_strain_nodal_rst() {
    RowMatrix rst(8, 3);
    rst << Precision(-1), Precision(-1), Precision(0),
           Precision( 1), Precision(-1), Precision(0),
           Precision( 1), Precision( 1), Precision(0),
           Precision(-1), Precision( 1), Precision(0),
           Precision( 0), Precision(-1), Precision(0),
           Precision( 1), Precision( 0), Precision(0),
           Precision( 0), Precision( 1), Precision(0),
           Precision(-1), Precision( 0), Precision(0);
    return rst;
}

FRTShellS8::VecN FRTShellS8::shape_function(Precision r, Precision s) const {
    VecN shape;

    // Corner shape functions
    shape(0) = Precision(0.25) * (Precision(1) - r) * (Precision(1) - s)
             * (-Precision(1) - r - s);
    shape(1) = Precision(0.25) * (Precision(1) + r) * (Precision(1) - s)
             * (-Precision(1) + r - s);
    shape(2) = Precision(0.25) * (Precision(1) + r) * (Precision(1) + s)
             * (-Precision(1) + r + s);
    shape(3) = Precision(0.25) * (Precision(1) - r) * (Precision(1) + s)
             * (-Precision(1) - r + s);

    // Midside shape functions
    shape(4) = Precision(0.5) * (Precision(1) - r * r) * (Precision(1) - s);
    shape(5) = Precision(0.5) * (Precision(1) + r) * (Precision(1) - s * s);
    shape(6) = Precision(0.5) * (Precision(1) - r * r) * (Precision(1) + s);
    shape(7) = Precision(0.5) * (Precision(1) - r) * (Precision(1) - s * s);

    return shape;
}

FRTShellS8::MatN2 FRTShellS8::shape_derivative(Precision r, Precision s) const {
    MatN2 derivative;

    // Derivatives of the four corner shape functions
    derivative(0, 0) = Precision(0.25) * (-Precision(2) * r - s) * (s - Precision(1));
    derivative(0, 1) = Precision(0.25) * (-r - Precision(2) * s) * (r - Precision(1));

    derivative(1, 0) = Precision(0.25) * (-Precision(2) * r + s) * (s - Precision(1));
    derivative(1, 1) = Precision(0.25) * (-r + Precision(2) * s) * (r + Precision(1));

    derivative(2, 0) = Precision(0.25) * (Precision(2) * r + s) * (s + Precision(1));
    derivative(2, 1) = Precision(0.25) * (r + Precision(1)) * (r + Precision(2) * s);

    derivative(3, 0) = Precision(0.25) * (Precision(2) * r - s) * (s + Precision(1));
    derivative(3, 1) = Precision(0.25) * (r - Precision(1)) * (r - Precision(2) * s);

    // Derivatives of the four midside shape functions
    derivative(4, 0) = r * (s - Precision(1));
    derivative(4, 1) = Precision(0.5) * (r * r - Precision(1));
    derivative(5, 0) = Precision(0.5) * (Precision(1) - s * s);
    derivative(5, 1) = -s * (r + Precision(1));
    derivative(6, 0) = -r * (Precision(1) + s);
    derivative(6, 1) = Precision(0.5) * (Precision(1) - r * r);
    derivative(7, 0) = Precision(0.5) * (s * s - Precision(1));
    derivative(7, 1) = s * (r - Precision(1));

    return derivative;
}

FRTShellS8::MatN2 FRTShellS8::node_coords_natural() const {
    MatN2 coordinates;
    coordinates << Precision(-1), Precision(-1),
                   Precision( 1), Precision(-1),
                   Precision( 1), Precision( 1),
                   Precision(-1), Precision( 1),
                   Precision( 0), Precision(-1),
                   Precision( 1), Precision( 0),
                   Precision( 0), Precision( 1),
                   Precision(-1), Precision( 0);
    return coordinates;
}

/**
 * Returns the classical MITC8 sampling positions.
 *
 * Points `0..7` are the eight in-layer interpolation points obtained by
 * scaling the serendipity nodal coordinates with `a = 1/sqrt(3)`. Points
 * `8..13` belong to the assumed `rt` shear field and points `14..19`
 * belong to the assumed `st` shear field. The fifth coefficient of each shear
 * interpolation is the mean of its two interior sampling values.
 */
std::vector<Vec2> FRTShellS8::tying_point_coordinates() const {
    const Precision a = Precision(1) / std::sqrt(Precision(3));

    return {
        // In-layer strain interpolation points
        Vec2(-a, -a),
        Vec2( a, -a),
        Vec2( a,  a),
        Vec2(-a,  a),
        Vec2(Precision(0), -a),
        Vec2( a, Precision(0)),
        Vec2(Precision(0),  a),
        Vec2(-a, Precision(0)),

        // rt shear boundary points and the RA/RB averaging points
        Vec2( a, Precision( 1)),
        Vec2(-a, Precision( 1)),
        Vec2(-a, Precision(-1)),
        Vec2( a, Precision(-1)),
        Vec2( a, Precision(0)),
        Vec2(-a, Precision(0)),

        // st shear boundary points and the SA/SB averaging points
        Vec2(Precision( 1),  a),
        Vec2(Precision( 1), -a),
        Vec2(Precision(-1), -a),
        Vec2(Precision(-1),  a),
        Vec2(Precision(0),  a),
        Vec2(Precision(0), -a)
    };
}

/**
 * Applies the classical MITC8 assumed in-layer and shear strain fields.
 *
 * The in-layer interpolation is performed on complete physical tangent
 * tensors. This is essential for curved or distorted elements because tensor
 * components sampled in different local bases cannot be averaged directly.
 * At the four midside points, the directly sampled edge-normal component is
 * combined with the adjacent corner tensors in an orthogonalized covariant
 * basis, following the original MITC8 construction.
 *
 * The `rt` and `st` transverse-shear fields are interpolated directly in
 * covariant natural components. Their fifth interpolation values are the means
 * at the RA/RB and SA/SB sampling pairs. The common shell kernel performs the
 * pointwise natural-to-local transformation after this function returns.
 */
void FRTShellS8::apply_mitc_natural(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat,
    Vec6NMat*             G_nat
) const {
    logging::error(data.tying_strain_nat.size() == 20,
                   "FRTShellS8: expected 20 MITC8 tying points");

    constexpr Index epsilon_start = 0;
    constexpr Index kappa_start   = 3;
    constexpr Index gamma_r3      = 6;
    constexpr Index gamma_s3      = 7;

    const auto& tying_points = reference_data().tying_points;
    const Precision a        = Precision(1) / std::sqrt(Precision(3));

    // Return the two contravariant surface vectors associated with the natural
    // coordinates of one reference point
    const auto dual_vectors = [](const ReferencePoint& sample) {
        std::array<Vec3, 2> dual;
        dual[0] = sample.invA(0, 0) * sample.e1 + sample.invA(1, 0) * sample.e2;
        dual[1] = sample.invA(0, 1) * sample.e1 + sample.invA(1, 1) * sample.e2;
        return dual;
    };

    // Map natural covariant in-plane components at a sampling point to the
    // engineering components measured along two supplied physical directions
    const auto natural_to_directions = [&](const ReferencePoint& sample,
                                           const Vec3&           direction_1,
                                           const Vec3&           direction_2) {
        const auto dual = dual_vectors(sample);

        const Precision u_r = direction_1.dot(dual[0]);
        const Precision u_s = direction_1.dot(dual[1]);
        const Precision v_r = direction_2.dot(dual[0]);
        const Precision v_s = direction_2.dot(dual[1]);

        StaticMatrix<3, 3> transformation;
        transformation << u_r*u_r,              u_s*u_s,              u_r*u_s,
                          v_r*v_r,              v_s*v_s,              v_r*v_s,
                          Precision(2)*u_r*v_r, Precision(2)*u_s*v_s, u_r*v_s + u_s*v_r;
        return transformation;
    };

    // Map covariant components in an arbitrary physical tangent basis to the
    // engineering components measured along two target directions
    const auto basis_to_directions = [](const Vec3& basis_1,
                                        const Vec3& basis_2,
                                        const Vec3& direction_1,
                                        const Vec3& direction_2) {
        Mat2 metric;
        metric << basis_1.dot(basis_1), basis_1.dot(basis_2),
                  basis_2.dot(basis_1), basis_2.dot(basis_2);

        const Precision determinant = metric.determinant();
        logging::error(std::abs(determinant) > Precision(1e-14),
                       "FRTShellS8: singular auxiliary MITC8 tangent basis");

        const Mat2 inverse = metric.inverse();
        const Vec3 dual_1  = inverse(0, 0) * basis_1 + inverse(1, 0) * basis_2;
        const Vec3 dual_2  = inverse(0, 1) * basis_1 + inverse(1, 1) * basis_2;

        const Precision u_1 = direction_1.dot(dual_1);
        const Precision u_2 = direction_1.dot(dual_2);
        const Precision v_1 = direction_2.dot(dual_1);
        const Precision v_2 = direction_2.dot(dual_2);

        StaticMatrix<3, 3> transformation;
        transformation << u_1*u_1,              u_2*u_2,              u_1*u_2,
                          v_1*v_1,              v_2*v_2,              v_1*v_2,
                          Precision(2)*u_1*v_1, Precision(2)*u_2*v_2, u_1*v_2 + u_2*v_1;
        return transformation;
    };

    // Transform three strain Hessians with one engineering tensor map
    const auto transform_G = [](const StaticMatrix<3, 3>& transformation,
                                const Vec6NMat&           source,
                                Index                     start) {
        std::array<Mat6N, 3> result;

        for (auto& value : result) {
            value.setZero();
        }

        for (Index target = 0; target < 3; ++target) {
            for (Index component = 0; component < 3; ++component) {
                result[target] += transformation(target, component) * source[start + component];
            }
        }

        return result;
    };

    // Construct the eight physical in-layer sampling tensors in the local
    // basis of the current integration point
    auto interpolate_in_layer_values = [&](Index start) {
        std::array<StaticVector<3>, 8> sampled;

        // Corner sampling tensors are taken directly from the compatible field
        for (Index corner = 0; corner < 4; ++corner) {
            const ReferencePoint& sample = tying_points[static_cast<std::size_t>(corner)];
            const StaticMatrix<3, 3> transformation =
                natural_to_directions(sample, point.e1, point.e2);

            sampled[corner] = transformation
                            * data.tying_strain_nat[static_cast<std::size_t>(corner)]
                                  .template segment<3>(start);
        }

        // Build one midside tensor. For bottom/top edges the directly evaluated
        // ss component is retained; for right/left edges the rr component is
        // retained. All remaining components follow from the adjacent corners.
        const auto build_midside = [&](Index midside,
                                       Index corner_1,
                                       Index corner_2,
                                       bool  keep_ss) {
            const ReferencePoint& middle = tying_points[static_cast<std::size_t>(midside)];

            Vec3 basis_1;
            Vec3 basis_2;

            if (keep_ss) {
                const Precision alpha = middle.X_r.dot(middle.X_s) / middle.X_s.squaredNorm();
                basis_1 = middle.X_r - alpha * middle.X_s;
                basis_2 = middle.X_s;
            } else {
                const Precision beta = middle.X_r.dot(middle.X_s) / middle.X_r.squaredNorm();
                basis_1 = middle.X_r;
                basis_2 = middle.X_s - beta * middle.X_r;
            }

            const StaticMatrix<3, 3> corner_1_to_basis =
                natural_to_directions(tying_points[static_cast<std::size_t>(corner_1)], basis_1, basis_2);
            const StaticMatrix<3, 3> corner_2_to_basis =
                natural_to_directions(tying_points[static_cast<std::size_t>(corner_2)], basis_1, basis_2);
            const StaticMatrix<3, 3> middle_to_basis =
                natural_to_directions(middle, basis_1, basis_2);

            const StaticVector<3> average = Precision(0.5)
                * (corner_1_to_basis
                   * data.tying_strain_nat[static_cast<std::size_t>(corner_1)].template segment<3>(start)
                 + corner_2_to_basis
                   * data.tying_strain_nat[static_cast<std::size_t>(corner_2)].template segment<3>(start));
            const StaticVector<3> direct = middle_to_basis
                * data.tying_strain_nat[static_cast<std::size_t>(midside)].template segment<3>(start);

            StaticVector<3> components;

            if (keep_ss) {
                components << average(0), direct(1), average(2);
            } else {
                components << direct(0), average(1), average(2);
            }

            const StaticMatrix<3, 3> to_target =
                basis_to_directions(basis_1, basis_2, point.e1, point.e2);
            const StaticVector<3> result = to_target * components;
            return result;
        };

        sampled[4] = build_midside(4, 0, 1, true);
        sampled[5] = build_midside(5, 1, 2, false);
        sampled[6] = build_midside(6, 3, 2, true);
        sampled[7] = build_midside(7, 0, 3, false);

        // Interpolate the eight physical tensors with serendipity functions on
        // the reduced square and convert the target-local components back to
        // natural covariant components
        const VecN weights = shape_function(point.r / a, point.s / a);
        StaticVector<3> local = StaticVector<3>::Zero();

        for (Index sample = 0; sample < 8; ++sample) {
            local += weights(sample) * sampled[sample];
        }

        const StaticMatrix<3, 3> natural_to_local =
            natural_to_directions(point, point.e1, point.e2);
        strain_nat.template segment<3>(start) = natural_to_local.inverse() * local;
    };

    interpolate_in_layer_values(epsilon_start);
    interpolate_in_layer_values(kappa_start);

    if (B_nat) {
        auto interpolate_in_layer_B = [&](Index start) {
            std::array<Mat3x6N, 8> sampled;

            for (Index corner = 0; corner < 4; ++corner) {
                const ReferencePoint& sample = tying_points[static_cast<std::size_t>(corner)];
                const StaticMatrix<3, 3> transformation =
                    natural_to_directions(sample, point.e1, point.e2);

                sampled[corner] = transformation
                    * data.tying_B_nat[static_cast<std::size_t>(corner)]
                          .template block<3, num_dofs>(start, 0);
            }

            const auto build_midside = [&](Index midside,
                                           Index corner_1,
                                           Index corner_2,
                                           bool  keep_ss) {
                const ReferencePoint& middle = tying_points[static_cast<std::size_t>(midside)];

                Vec3 basis_1;
                Vec3 basis_2;

                if (keep_ss) {
                    const Precision alpha = middle.X_r.dot(middle.X_s) / middle.X_s.squaredNorm();
                    basis_1 = middle.X_r - alpha * middle.X_s;
                    basis_2 = middle.X_s;
                } else {
                    const Precision beta = middle.X_r.dot(middle.X_s) / middle.X_r.squaredNorm();
                    basis_1 = middle.X_r;
                    basis_2 = middle.X_s - beta * middle.X_r;
                }

                const StaticMatrix<3, 3> corner_1_to_basis =
                    natural_to_directions(tying_points[static_cast<std::size_t>(corner_1)], basis_1, basis_2);
                const StaticMatrix<3, 3> corner_2_to_basis =
                    natural_to_directions(tying_points[static_cast<std::size_t>(corner_2)], basis_1, basis_2);
                const StaticMatrix<3, 3> middle_to_basis =
                    natural_to_directions(middle, basis_1, basis_2);

                const Mat3x6N average = Precision(0.5)
                    * (corner_1_to_basis
                       * data.tying_B_nat[static_cast<std::size_t>(corner_1)]
                             .template block<3, num_dofs>(start, 0)
                     + corner_2_to_basis
                       * data.tying_B_nat[static_cast<std::size_t>(corner_2)]
                             .template block<3, num_dofs>(start, 0));
                const Mat3x6N direct = middle_to_basis
                    * data.tying_B_nat[static_cast<std::size_t>(midside)]
                          .template block<3, num_dofs>(start, 0);

                Mat3x6N components;
                components.row(0) = keep_ss ? average.row(0) : direct.row(0);
                components.row(1) = keep_ss ? direct.row(1)  : average.row(1);
                components.row(2) = average.row(2);

                const StaticMatrix<3, 3> to_target =
                    basis_to_directions(basis_1, basis_2, point.e1, point.e2);
                const Mat3x6N result = to_target * components;
                return result;
            };

            sampled[4] = build_midside(4, 0, 1, true);
            sampled[5] = build_midside(5, 1, 2, false);
            sampled[6] = build_midside(6, 3, 2, true);
            sampled[7] = build_midside(7, 0, 3, false);

            const VecN weights = shape_function(point.r / a, point.s / a);
            Mat3x6N local = Mat3x6N::Zero();

            for (Index sample = 0; sample < 8; ++sample) {
                local += weights(sample) * sampled[sample];
            }

            const StaticMatrix<3, 3> natural_to_local =
                natural_to_directions(point, point.e1, point.e2);
            B_nat->template block<3, num_dofs>(start, 0) = natural_to_local.inverse() * local;
        };

        interpolate_in_layer_B(epsilon_start);
        interpolate_in_layer_B(kappa_start);
    }

    if (G_nat) {
        auto interpolate_in_layer_G = [&](Index start) {
            std::array<std::array<Mat6N, 3>, 8> sampled;

            for (Index corner = 0; corner < 4; ++corner) {
                const ReferencePoint& sample = tying_points[static_cast<std::size_t>(corner)];
                const StaticMatrix<3, 3> transformation =
                    natural_to_directions(sample, point.e1, point.e2);

                sampled[corner] = transform_G(
                    transformation,
                    data.tying_G_nat[static_cast<std::size_t>(corner)],
                    start
                );
            }

            const auto build_midside = [&](Index midside,
                                           Index corner_1,
                                           Index corner_2,
                                           bool  keep_ss) {
                const ReferencePoint& middle = tying_points[static_cast<std::size_t>(midside)];

                Vec3 basis_1;
                Vec3 basis_2;

                if (keep_ss) {
                    const Precision alpha = middle.X_r.dot(middle.X_s) / middle.X_s.squaredNorm();
                    basis_1 = middle.X_r - alpha * middle.X_s;
                    basis_2 = middle.X_s;
                } else {
                    const Precision beta = middle.X_r.dot(middle.X_s) / middle.X_r.squaredNorm();
                    basis_1 = middle.X_r;
                    basis_2 = middle.X_s - beta * middle.X_r;
                }

                const StaticMatrix<3, 3> corner_1_to_basis =
                    natural_to_directions(tying_points[static_cast<std::size_t>(corner_1)], basis_1, basis_2);
                const StaticMatrix<3, 3> corner_2_to_basis =
                    natural_to_directions(tying_points[static_cast<std::size_t>(corner_2)], basis_1, basis_2);
                const StaticMatrix<3, 3> middle_to_basis =
                    natural_to_directions(middle, basis_1, basis_2);

                const auto corner_1_G = transform_G(
                    corner_1_to_basis,
                    data.tying_G_nat[static_cast<std::size_t>(corner_1)],
                    start
                );
                const auto corner_2_G = transform_G(
                    corner_2_to_basis,
                    data.tying_G_nat[static_cast<std::size_t>(corner_2)],
                    start
                );
                const auto middle_G = transform_G(
                    middle_to_basis,
                    data.tying_G_nat[static_cast<std::size_t>(midside)],
                    start
                );

                std::array<Mat6N, 3> components;
                components[0] = keep_ss
                    ? Precision(0.5) * (corner_1_G[0] + corner_2_G[0])
                    : middle_G[0];
                components[1] = keep_ss
                    ? middle_G[1]
                    : Precision(0.5) * (corner_1_G[1] + corner_2_G[1]);
                components[2] = Precision(0.5) * (corner_1_G[2] + corner_2_G[2]);

                const StaticMatrix<3, 3> to_target =
                    basis_to_directions(basis_1, basis_2, point.e1, point.e2);
                std::array<Mat6N, 3> result;

                for (auto& value : result) {
                    value.setZero();
                }

                for (Index target = 0; target < 3; ++target) {
                    for (Index component = 0; component < 3; ++component) {
                        result[target] += to_target(target, component) * components[component];
                    }
                }

                return result;
            };

            sampled[4] = build_midside(4, 0, 1, true);
            sampled[5] = build_midside(5, 1, 2, false);
            sampled[6] = build_midside(6, 3, 2, true);
            sampled[7] = build_midside(7, 0, 3, false);

            const VecN weights = shape_function(point.r / a, point.s / a);
            std::array<Mat6N, 3> local;

            for (auto& value : local) {
                value.setZero();
            }

            for (Index sample = 0; sample < 8; ++sample) {
                for (Index component = 0; component < 3; ++component) {
                    local[component] += weights(sample) * sampled[sample][component];
                }
            }

            const StaticMatrix<3, 3> local_to_natural =
                natural_to_directions(point, point.e1, point.e2).inverse();

            for (Index target = 0; target < 3; ++target) {
                (*G_nat)[start + target].setZero();

                for (Index component = 0; component < 3; ++component) {
                    (*G_nat)[start + target] +=
                        local_to_natural(target, component) * local[component];
                }
            }
        };

        interpolate_in_layer_G(epsilon_start);
        interpolate_in_layer_G(kappa_start);
    }

    // Interpolate the covariant rt shear component. The four boundary
    // functions and the fifth internal function follow the classical MITC8
    // construction. The fifth coefficient is the mean of the RA/RB values.
    const Precision h_rt_5 = Precision(1) - point.s * point.s;
    const Precision h_rt_1 = Precision(0.25) * (Precision(1) + point.r / a)
                           * (Precision(1) + point.s) - Precision(0.25) * h_rt_5;
    const Precision h_rt_2 = Precision(0.25) * (Precision(1) - point.r / a)
                           * (Precision(1) + point.s) - Precision(0.25) * h_rt_5;
    const Precision h_rt_3 = Precision(0.25) * (Precision(1) - point.r / a)
                           * (Precision(1) - point.s) - Precision(0.25) * h_rt_5;
    const Precision h_rt_4 = Precision(0.25) * (Precision(1) + point.r / a)
                           * (Precision(1) - point.s) - Precision(0.25) * h_rt_5;

    const Precision rt_mean = Precision(0.5)
                            * (data.tying_strain_nat[12](gamma_r3)
                             + data.tying_strain_nat[13](gamma_r3));

    strain_nat(gamma_r3) =
          h_rt_1 * data.tying_strain_nat[ 8](gamma_r3)
        + h_rt_2 * data.tying_strain_nat[ 9](gamma_r3)
        + h_rt_3 * data.tying_strain_nat[10](gamma_r3)
        + h_rt_4 * data.tying_strain_nat[11](gamma_r3)
        + h_rt_5 * rt_mean;

    // Interpolate the covariant st shear component. The fifth coefficient is
    // the mean of the SA/SB values.
    const Precision h_st_5 = Precision(1) - point.r * point.r;
    const Precision h_st_1 = Precision(0.25) * (Precision(1) + point.s / a)
                           * (Precision(1) + point.r) - Precision(0.25) * h_st_5;
    const Precision h_st_2 = Precision(0.25) * (Precision(1) + point.s / a)
                           * (Precision(1) - point.r) - Precision(0.25) * h_st_5;
    const Precision h_st_3 = Precision(0.25) * (Precision(1) - point.s / a)
                           * (Precision(1) - point.r) - Precision(0.25) * h_st_5;
    const Precision h_st_4 = Precision(0.25) * (Precision(1) - point.s / a)
                           * (Precision(1) + point.r) - Precision(0.25) * h_st_5;

    const Precision st_mean = Precision(0.5)
                            * (data.tying_strain_nat[18](gamma_s3)
                             + data.tying_strain_nat[19](gamma_s3));

    strain_nat(gamma_s3) =
          h_st_1 * data.tying_strain_nat[14](gamma_s3)
        + h_st_2 * data.tying_strain_nat[17](gamma_s3)
        + h_st_3 * data.tying_strain_nat[16](gamma_s3)
        + h_st_4 * data.tying_strain_nat[15](gamma_s3)
        + h_st_5 * st_mean;

    if (B_nat) {
        const StaticMatrix<1, num_dofs> rt_B_mean =
            Precision(0.5) * (data.tying_B_nat[12].row(gamma_r3)
                            + data.tying_B_nat[13].row(gamma_r3));
        const StaticMatrix<1, num_dofs> st_B_mean =
            Precision(0.5) * (data.tying_B_nat[18].row(gamma_s3)
                            + data.tying_B_nat[19].row(gamma_s3));

        B_nat->row(gamma_r3) =
              h_rt_1 * data.tying_B_nat[ 8].row(gamma_r3)
            + h_rt_2 * data.tying_B_nat[ 9].row(gamma_r3)
            + h_rt_3 * data.tying_B_nat[10].row(gamma_r3)
            + h_rt_4 * data.tying_B_nat[11].row(gamma_r3)
            + h_rt_5 * rt_B_mean;

        B_nat->row(gamma_s3) =
              h_st_1 * data.tying_B_nat[14].row(gamma_s3)
            + h_st_2 * data.tying_B_nat[17].row(gamma_s3)
            + h_st_3 * data.tying_B_nat[16].row(gamma_s3)
            + h_st_4 * data.tying_B_nat[15].row(gamma_s3)
            + h_st_5 * st_B_mean;
    }

    if (G_nat) {
        const Mat6N rt_G_mean = Precision(0.5)
                              * (data.tying_G_nat[12][gamma_r3]
                               + data.tying_G_nat[13][gamma_r3]);
        const Mat6N st_G_mean = Precision(0.5)
                              * (data.tying_G_nat[18][gamma_s3]
                               + data.tying_G_nat[19][gamma_s3]);

        (*G_nat)[gamma_r3] =
              h_rt_1 * data.tying_G_nat[ 8][gamma_r3]
            + h_rt_2 * data.tying_G_nat[ 9][gamma_r3]
            + h_rt_3 * data.tying_G_nat[10][gamma_r3]
            + h_rt_4 * data.tying_G_nat[11][gamma_r3]
            + h_rt_5 * rt_G_mean;

        (*G_nat)[gamma_s3] =
              h_st_1 * data.tying_G_nat[14][gamma_s3]
            + h_st_2 * data.tying_G_nat[17][gamma_s3]
            + h_st_3 * data.tying_G_nat[16][gamma_s3]
            + h_st_4 * data.tying_G_nat[15][gamma_s3]
            + h_st_5 * st_G_mean;
    }

}

} // namespace fem::model
