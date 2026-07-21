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
 * applied identically to the strain values and B rows. Geometric-tangent
 * weights use the corresponding transposed interpolation.
 *
 * @see FRTShellS8
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell_s8.h"

#include "../../core/logging.h"

#include <cmath>

namespace fem::model {

/**
 * Constructs the eight-node MITC8 finite-rotation shell.
 *
 * @param id Element identifier.
 * @param nodes Eight nodal identifiers in quadratic serendipity ordering.
 */
FRTShellS8::FRTShellS8(ID id, const std::array<ID, 8>& nodes)
    : FRTShell<8>       (id, nodes),
      integration_scheme_(math::quadrature::Domain::DOMAIN_ISO_QUAD,
                          math::quadrature::Order::ORDER_QUINTIC) {}

/**
 * Returns the public element type identifier.
 *
 * @return String identifier `MITC8FRT`.
 */
std::string FRTShellS8::type_name() const {
    return "MITC8FRT";
}

/**
 * Returns an oriented quadratic serendipity surface representation.
 *
 * The reverse side changes both corner and midside ordering consistently.
 *
 * @param surface_id Requested positive or negative shell side.
 * @return Oriented eight-node surface object.
 */
std::shared_ptr<SurfaceInterface> FRTShellS8::surface(int surface_id) {
    return std::make_shared<Surface8>(
        surface_id == 1
            ? std::array<ID, 8>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3],
                                this->nodes()[4], this->nodes()[5], this->nodes()[6], this->nodes()[7]}
            : std::array<ID, 8>{this->nodes()[0], this->nodes()[3], this->nodes()[2], this->nodes()[1],
                                this->nodes()[7], this->nodes()[6], this->nodes()[5], this->nodes()[4]}
    );
}

/**
 * Returns the full quadrilateral area quadrature used by the element.
 *
 * @return Quintic rule corresponding to full `3 x 3` integration.
 */
const math::quadrature::Quadrature& FRTShellS8::integration_scheme() const {
    return integration_scheme_;
}

/**
 * Returns all numerical integration-point coordinates for result storage.
 *
 * @return Matrix of `(r,s,t)` integration-point coordinates.
 */
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

/**
 * Returns the eight natural nodal coordinates for result recovery.
 *
 * @return Matrix of corner and midside coordinates in `(r,s,t)` format.
 */
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

/**
 * Evaluates the eight quadratic serendipity shape functions.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Eight serendipity shape-function values.
 */
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

/**
 * Evaluates the first natural derivatives of the serendipity interpolation.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Matrix whose rows contain `[dN_i/dr,dN_i/ds]`.
 */
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

/**
 * Returns the natural coordinates of the four corner and four midside nodes.
 *
 * @return Eight serendipity node coordinates on `[-1,1]^2`.
 */
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
 *
 * @return Ordered natural tying-point coordinates.
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
 *
 * @param data Active tying-point strain and B data.
 * @param point Target natural shell point.
 * @param strain_nat Compatible natural strain vector to modify.
 * @param B_nat Optional compatible natural B matrix to modify.
 */
void FRTShellS8::apply_mitc_natural(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat
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
        dual[0] = sample.invJ(0, 0) * sample.basis.col(0) + sample.invJ(1, 0) * sample.basis.col(1);
        dual[1] = sample.invJ(0, 1) * sample.basis.col(0) + sample.invJ(1, 1) * sample.basis.col(1);
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

    // Construct the eight physical in-layer sampling tensors in the local
    // basis of the current integration point
    auto interpolate_in_layer_values = [&](Index start) {
        std::array<StaticVector<3>, 8> sampled;

        // Corner sampling tensors are taken directly from the compatible field
        for (Index corner = 0; corner < 4; ++corner) {
            const ReferencePoint& sample = tying_points[static_cast<std::size_t>(corner)];
            const StaticMatrix<3, 3> transformation =
                natural_to_directions(sample, point.basis.col(0), point.basis.col(1));

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
                const Precision alpha = middle.X_rs.col(0).dot(middle.X_rs.col(1)) / middle.X_rs.col(1).squaredNorm();
                basis_1 = middle.X_rs.col(0) - alpha * middle.X_rs.col(1);
                basis_2 = middle.X_rs.col(1);
            } else {
                const Precision beta = middle.X_rs.col(0).dot(middle.X_rs.col(1)) / middle.X_rs.col(0).squaredNorm();
                basis_1 = middle.X_rs.col(0);
                basis_2 = middle.X_rs.col(1) - beta * middle.X_rs.col(0);
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
                basis_to_directions(basis_1, basis_2, point.basis.col(0), point.basis.col(1));
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
            natural_to_directions(point, point.basis.col(0), point.basis.col(1));
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
                    natural_to_directions(sample, point.basis.col(0), point.basis.col(1));

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
                    const Precision alpha = middle.X_rs.col(0).dot(middle.X_rs.col(1)) / middle.X_rs.col(1).squaredNorm();
                    basis_1 = middle.X_rs.col(0) - alpha * middle.X_rs.col(1);
                    basis_2 = middle.X_rs.col(1);
                } else {
                    const Precision beta = middle.X_rs.col(0).dot(middle.X_rs.col(1)) / middle.X_rs.col(0).squaredNorm();
                    basis_1 = middle.X_rs.col(0);
                    basis_2 = middle.X_rs.col(1) - beta * middle.X_rs.col(0);
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
                    basis_to_directions(basis_1, basis_2, point.basis.col(0), point.basis.col(1));
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
                natural_to_directions(point, point.basis.col(0), point.basis.col(1));
            B_nat->template block<3, num_dofs>(start, 0) = natural_to_local.inverse() * local;
        };

        interpolate_in_layer_B(epsilon_start);
        interpolate_in_layer_B(kappa_start);
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

}

/**
 * Applies the exact transpose of the complete classical MITC8 interpolation.
 *
 * In-layer tensor weights are transformed between the physical sampling bases
 * and distributed to the eight compatible sampling tensors. The two assumed
 * shear weights are distributed to their boundary and interior averaging
 * samples. The caller supplies a zeroed thread-local span, so no dynamic
 * allocation occurs in the geometric-tangent loop.
 *
 * @param point Integration point defining all MITC8 interpolation functions.
 * @param assumed_weights Generalized natural weights after basis pull-back.
 * @param compatible_weights Compatible integration-point weights to increment.
 * @param tying_weights Zeroed compatible tying-point weights to increment.
 */
void FRTShellS8::pull_back_mitc_resultants(
    const ReferencePoint& point,
    const Vec8&           assumed_weights,
    Vec8&                 compatible_weights,
    Span<Vec8>            tying_weights
) const {
    (void) compatible_weights;

    constexpr Index epsilon_start = 0;
    constexpr Index kappa_start   = 3;
    constexpr Index gamma_r3      = 6;
    constexpr Index gamma_s3      = 7;

    const auto& tying_points = reference_data().tying_points;
    const Precision a        = Precision(1) / std::sqrt(Precision(3));

    const auto dual_vectors = [](const ReferencePoint& sample) {
        std::array<Vec3, 2> dual;
        dual[0] = sample.invJ(0, 0) * sample.basis.col(0) + sample.invJ(1, 0) * sample.basis.col(1);
        dual[1] = sample.invJ(0, 1) * sample.basis.col(0) + sample.invJ(1, 1) * sample.basis.col(1);
        return dual;
    };

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

    const auto pull_back_in_layer = [&](Index start) {
        std::array<StaticVector<3>, 8> sampled_weights;

        for (auto& value : sampled_weights) {
            value.setZero();
        }

        const StaticMatrix<3, 3> local_to_natural =
            natural_to_directions(point, point.basis.col(0), point.basis.col(1)).inverse();
        const StaticVector<3> local_weight =
            local_to_natural.transpose() * assumed_weights.template segment<3>(start);
        const VecN weights = shape_function(point.r / a, point.s / a);

        for (Index sample = 0; sample < 8; ++sample) {
            sampled_weights[sample] += weights(sample) * local_weight;
        }

        for (Index corner = 0; corner < 4; ++corner) {
            const ReferencePoint& sample = tying_points[static_cast<std::size_t>(corner)];
            const StaticMatrix<3, 3> transformation =
                natural_to_directions(sample, point.basis.col(0), point.basis.col(1));

            tying_weights[static_cast<std::size_t>(corner)].template segment<3>(start) +=
                transformation.transpose() * sampled_weights[corner];
        }

        const auto pull_back_midside = [&](Index midside,
                                           Index corner_1,
                                           Index corner_2,
                                           bool  keep_ss) {
            const ReferencePoint& middle = tying_points[static_cast<std::size_t>(midside)];

            Vec3 basis_1;
            Vec3 basis_2;

            if (keep_ss) {
                const Precision alpha = middle.X_rs.col(0).dot(middle.X_rs.col(1)) / middle.X_rs.col(1).squaredNorm();
                basis_1 = middle.X_rs.col(0) - alpha * middle.X_rs.col(1);
                basis_2 = middle.X_rs.col(1);
            } else {
                const Precision beta = middle.X_rs.col(0).dot(middle.X_rs.col(1)) / middle.X_rs.col(0).squaredNorm();
                basis_1 = middle.X_rs.col(0);
                basis_2 = middle.X_rs.col(1) - beta * middle.X_rs.col(0);
            }

            const StaticMatrix<3, 3> corner_1_to_basis =
                natural_to_directions(tying_points[static_cast<std::size_t>(corner_1)], basis_1, basis_2);
            const StaticMatrix<3, 3> corner_2_to_basis =
                natural_to_directions(tying_points[static_cast<std::size_t>(corner_2)], basis_1, basis_2);
            const StaticMatrix<3, 3> middle_to_basis =
                natural_to_directions(middle, basis_1, basis_2);
            const StaticMatrix<3, 3> to_target =
                basis_to_directions(basis_1, basis_2, point.basis.col(0), point.basis.col(1));

            const StaticVector<3> component_weights =
                to_target.transpose() * sampled_weights[midside];
            StaticVector<3> average_weights = StaticVector<3>::Zero();
            StaticVector<3> direct_weights  = StaticVector<3>::Zero();

            if (keep_ss) {
                average_weights(0) += component_weights(0);
                direct_weights (1) += component_weights(1);
                average_weights(2) += component_weights(2);
            } else {
                direct_weights (0) += component_weights(0);
                average_weights(1) += component_weights(1);
                average_weights(2) += component_weights(2);
            }

            tying_weights[static_cast<std::size_t>(corner_1)].template segment<3>(start) +=
                Precision(0.5) * corner_1_to_basis.transpose() * average_weights;
            tying_weights[static_cast<std::size_t>(corner_2)].template segment<3>(start) +=
                Precision(0.5) * corner_2_to_basis.transpose() * average_weights;
            tying_weights[static_cast<std::size_t>(midside)].template segment<3>(start) +=
                middle_to_basis.transpose() * direct_weights;
        };

        pull_back_midside(4, 0, 1, true);
        pull_back_midside(5, 1, 2, false);
        pull_back_midside(6, 3, 2, true);
        pull_back_midside(7, 0, 3, false);
    };

    pull_back_in_layer(epsilon_start);
    pull_back_in_layer(kappa_start);

    const Precision h_rt_5 = Precision(1) - point.s * point.s;
    const Precision h_rt_1 = Precision(0.25) * (Precision(1) + point.r / a)
                           * (Precision(1) + point.s) - Precision(0.25) * h_rt_5;
    const Precision h_rt_2 = Precision(0.25) * (Precision(1) - point.r / a)
                           * (Precision(1) + point.s) - Precision(0.25) * h_rt_5;
    const Precision h_rt_3 = Precision(0.25) * (Precision(1) - point.r / a)
                           * (Precision(1) - point.s) - Precision(0.25) * h_rt_5;
    const Precision h_rt_4 = Precision(0.25) * (Precision(1) + point.r / a)
                           * (Precision(1) - point.s) - Precision(0.25) * h_rt_5;

    tying_weights[ 8](gamma_r3) += h_rt_1 * assumed_weights(gamma_r3);
    tying_weights[ 9](gamma_r3) += h_rt_2 * assumed_weights(gamma_r3);
    tying_weights[10](gamma_r3) += h_rt_3 * assumed_weights(gamma_r3);
    tying_weights[11](gamma_r3) += h_rt_4 * assumed_weights(gamma_r3);
    tying_weights[12](gamma_r3) += Precision(0.5) * h_rt_5 * assumed_weights(gamma_r3);
    tying_weights[13](gamma_r3) += Precision(0.5) * h_rt_5 * assumed_weights(gamma_r3);

    const Precision h_st_5 = Precision(1) - point.r * point.r;
    const Precision h_st_1 = Precision(0.25) * (Precision(1) + point.s / a)
                           * (Precision(1) + point.r) - Precision(0.25) * h_st_5;
    const Precision h_st_2 = Precision(0.25) * (Precision(1) + point.s / a)
                           * (Precision(1) - point.r) - Precision(0.25) * h_st_5;
    const Precision h_st_3 = Precision(0.25) * (Precision(1) - point.s / a)
                           * (Precision(1) - point.r) - Precision(0.25) * h_st_5;
    const Precision h_st_4 = Precision(0.25) * (Precision(1) - point.s / a)
                           * (Precision(1) + point.r) - Precision(0.25) * h_st_5;

    tying_weights[14](gamma_s3) += h_st_1 * assumed_weights(gamma_s3);
    tying_weights[17](gamma_s3) += h_st_2 * assumed_weights(gamma_s3);
    tying_weights[16](gamma_s3) += h_st_3 * assumed_weights(gamma_s3);
    tying_weights[15](gamma_s3) += h_st_4 * assumed_weights(gamma_s3);
    tying_weights[18](gamma_s3) += Precision(0.5) * h_st_5 * assumed_weights(gamma_s3);
    tying_weights[19](gamma_s3) += Precision(0.5) * h_st_5 * assumed_weights(gamma_s3);
}

} // namespace fem::model
