#include "so3.h"

#include "vec_util.h"

#include <cmath>

namespace fem::so3 {

namespace {

using math::skew;

struct RotationCoefficients {
    Precision a     = Precision(0);
    Precision b     = Precision(0);
    Precision a_s   = Precision(0);
    Precision b_s   = Precision(0);
    Precision a_ss  = Precision(0);
    Precision b_ss  = Precision(0);
};


RotationCoefficients rotation_coefficients(Precision angle_squared) {
    RotationCoefficients coefficients;

    if (angle_squared < Precision(1e-4)) {
        const Precision s2 = angle_squared * angle_squared;
        const Precision s3 = s2 * angle_squared;

        coefficients.a    = Precision(1)
                          - angle_squared / Precision(6)
                          + s2            / Precision(120)
                          - s3            / Precision(5040);
        coefficients.b    = Precision(0.5)
                          - angle_squared / Precision(24)
                          + s2            / Precision(720)
                          - s3            / Precision(40320);
        coefficients.a_s  = -Precision(1) / Precision(6)
                          + angle_squared / Precision(60)
                          - s2            / Precision(1680)
                          + s3            / Precision(90720);
        coefficients.b_s  = -Precision(1) / Precision(24)
                          + angle_squared / Precision(360)
                          - s2            / Precision(13440)
                          + s3            / Precision(907200);
        coefficients.a_ss = Precision(1) / Precision(60)
                          - angle_squared / Precision(840)
                          + s2            / Precision(30240);
        coefficients.b_ss = Precision(1) / Precision(360)
                          - angle_squared / Precision(6720)
                          + s2            / Precision(302400);

        return coefficients;
    }

    const Precision angle        = std::sqrt(angle_squared);
    const Precision angle_cubed  = angle_squared * angle;
    const Precision angle_fourth = angle_squared * angle_squared;
    const Precision angle_fifth  = angle_fourth * angle;
    const Precision angle_sixth  = angle_fourth * angle_squared;
    const Precision sin_angle    = std::sin(angle);
    const Precision cos_angle    = std::cos(angle);

    coefficients.a    = sin_angle / angle;
    coefficients.b    = (Precision(1) - cos_angle) / angle_squared;
    coefficients.a_s  = (angle * cos_angle - sin_angle)
                      / (Precision(2) * angle_cubed);
    coefficients.b_s  = (Precision(0.5) * angle * sin_angle
                       + cos_angle - Precision(1))
                      / angle_fourth;
    coefficients.a_ss = (-angle_squared * sin_angle
                       - Precision(3) * angle * cos_angle
                       + Precision(3) * sin_angle)
                      / (Precision(4) * angle_fifth);
    coefficients.b_ss = (angle_squared * cos_angle
                       - Precision(5) * angle * sin_angle
                       - Precision(8) * cos_angle
                       + Precision(8))
                      / (Precision(4) * angle_sixth);

    return coefficients;
}


std::array<Mat3, 3> rotation_generators() {
    return {
        skew(Vec3::UnitX()),
        skew(Vec3::UnitY()),
        skew(Vec3::UnitZ())
    };
}

} // namespace


Mat3 rotation_matrix(const Vec3& rotation_vector) {
    const auto coefficients = rotation_coefficients(rotation_vector.squaredNorm());

    const Mat3 omega         = skew(rotation_vector);
    const Mat3 omega_squared = omega * omega;

    return Mat3::Identity()
         + coefficients.a * omega
         + coefficients.b * omega_squared;
}


void rotation_matrix_first_derivatives(
    const Vec3&          rotation_vector,
    Mat3&                rotation,
    std::array<Mat3, 3>& first_derivatives
) {
    const Precision angle_squared = rotation_vector.squaredNorm();
    const auto      coefficients  = rotation_coefficients(angle_squared);
    const Mat3      omega         = skew(rotation_vector);
    const Mat3      omega_squared = omega * omega;
    const auto      generators    = rotation_generators();

    rotation = Mat3::Identity()
             + coefficients.a * omega
             + coefficients.b * omega_squared;

    for (Index i = 0; i < 3; ++i) {
        const Precision a_i = Precision(2)
                            * coefficients.a_s
                            * rotation_vector(i);
        const Precision b_i = Precision(2)
                            * coefficients.b_s
                            * rotation_vector(i);

        const Mat3 omega_squared_derivative = generators[i] * omega
                                            + omega * generators[i];

        first_derivatives[i] = a_i * omega
                             + coefficients.a * generators[i]
                             + b_i * omega_squared
                             + coefficients.b * omega_squared_derivative;
    }
}


void rotation_matrix_second_derivatives(
    const Vec3&                         rotation_vector,
    Mat3&                               rotation,
    std::array<Mat3, 3>&                first_derivatives,
    std::array<std::array<Mat3, 3>, 3>& second_derivatives
) {
    const Precision angle_squared = rotation_vector.squaredNorm();
    const auto      coefficients  = rotation_coefficients(angle_squared);
    const Mat3      omega         = skew(rotation_vector);
    const Mat3      omega_squared = omega * omega;
    const auto      generators    = rotation_generators();

    std::array<Mat3, 3> omega_squared_derivatives;

    rotation = Mat3::Identity()
             + coefficients.a * omega
             + coefficients.b * omega_squared;

    for (Index i = 0; i < 3; ++i) {
        const Precision a_i = Precision(2)
                            * coefficients.a_s
                            * rotation_vector(i);
        const Precision b_i = Precision(2)
                            * coefficients.b_s
                            * rotation_vector(i);

        omega_squared_derivatives[i] = generators[i] * omega
                                     + omega * generators[i];

        first_derivatives[i] = a_i * omega
                             + coefficients.a * generators[i]
                             + b_i * omega_squared
                             + coefficients.b * omega_squared_derivatives[i];
    }

    for (Index i = 0; i < 3; ++i) {
        const Precision a_i = Precision(2)
                            * coefficients.a_s
                            * rotation_vector(i);
        const Precision b_i = Precision(2)
                            * coefficients.b_s
                            * rotation_vector(i);

        for (Index j = 0; j < 3; ++j) {
            const Precision delta = i == j ? Precision(1) : Precision(0);

            const Precision a_j = Precision(2)
                                * coefficients.a_s
                                * rotation_vector(j);
            const Precision b_j = Precision(2)
                                * coefficients.b_s
                                * rotation_vector(j);

            const Precision a_ij = Precision(4)
                                 * coefficients.a_ss
                                 * rotation_vector(i)
                                 * rotation_vector(j)
                                 + Precision(2)
                                 * coefficients.a_s
                                 * delta;
            const Precision b_ij = Precision(4)
                                 * coefficients.b_ss
                                 * rotation_vector(i)
                                 * rotation_vector(j)
                                 + Precision(2)
                                 * coefficients.b_s
                                 * delta;

            second_derivatives[i][j] =
                  a_ij * omega
                + a_i  * generators[j]
                + a_j  * generators[i]
                + b_ij * omega_squared
                + b_i  * omega_squared_derivatives[j]
                + b_j  * omega_squared_derivatives[i]
                + coefficients.b
                * (generators[i] * generators[j]
                 + generators[j] * generators[i]);
        }
    }
}

} // namespace fem::so3
