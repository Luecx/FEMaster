/**
 * @file stress.cpp
 * @brief Implements stress rotation utilities in Voigt notation.
 *
 * @see src/material/stress.h
 * @see src/cos/coordinate_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "stress.h"

#include <Eigen/LU>

namespace fem {

Stress::Stress(const Vec6& voigt)
    : voigt_(voigt) {}

Precision& Stress::operator[](Component component) {
    return voigt_(static_cast<int>(component));
}

Precision Stress::operator[](Component component) const {
    return voigt_(static_cast<int>(component));
}

const Vec6& Stress::voigt() const {
    return voigt_;
}

Vec6& Stress::voigt() {
    return voigt_;
}

Mat3 Stress::tensor() const {
    Mat3 tensor;
    tensor << voigt_(0), voigt_(5), voigt_(4),
              voigt_(5), voigt_(1), voigt_(3),
              voigt_(4), voigt_(3), voigt_(2);
    return tensor;
}

Stress Stress::transformed(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) const {
    return Stress{Stress::get_transformation_matrix(from_basis, to_basis) * voigt_};
}

Stress Stress::transform(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) const {
    return transformed(from_basis, to_basis);
}

Stress Stress::from_tensor(const Mat3& tensor) {
    Vec6 voigt;
    voigt << tensor(0, 0),
             tensor(1, 1),
             tensor(2, 2),
             tensor(1, 2),
             tensor(2, 0),
             tensor(0, 1);
    return Stress{voigt};
}

/**
 * @copydoc Stress::get_transformation_matrix
 */
Mat6 Stress::get_transformation_matrix(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) {
    /*
     * Basis convention:
     *
     * Each basis is assumed to store its basis vectors as columns expressed in
     * the same parent/global coordinate system.
     *
     * For a passive component transformation from from_basis to to_basis:
     *
     *     sigma_to = R * sigma_from * R^T
     *
     * with
     *
     *     R = to_basis^T * from_basis
     *
     * where R(i,j) = e_i_to dot e_j_from.
     */
    const Mat3 R = to_basis.transpose() * from_basis;

    const Precision R11 = R(0, 0);
    const Precision R12 = R(0, 1);
    const Precision R13 = R(0, 2);

    const Precision R21 = R(1, 0);
    const Precision R22 = R(1, 1);
    const Precision R23 = R(1, 2);

    const Precision R31 = R(2, 0);
    const Precision R32 = R(2, 1);
    const Precision R33 = R(2, 2);

    StaticMatrix<6, 6> T_sigma;

    T_sigma <<
        R11 * R11, R12 * R12, R13 * R13, 2 * R12 * R13, 2 * R11 * R13, 2 * R11 * R12,
        R21 * R21, R22 * R22, R23 * R23, 2 * R22 * R23, 2 * R21 * R23, 2 * R21 * R22,
        R31 * R31, R32 * R32, R33 * R33, 2 * R32 * R33, 2 * R31 * R33, 2 * R31 * R32,
        R21 * R31, R22 * R32, R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
        R11 * R31, R12 * R32, R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
        R11 * R21, R12 * R22, R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;

    return T_sigma;
}

CauchyStress CauchyStress::transformed(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) const {
    return CauchyStress{Stress::get_transformation_matrix(from_basis, to_basis) * voigt_};
}

CauchyStress CauchyStress::from_tensor(const Mat3& tensor) {
    return CauchyStress{Stress::from_tensor(tensor).voigt()};
}

PK2Stress PK2Stress::transformed(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) const {
    return PK2Stress{Stress::get_transformation_matrix(from_basis, to_basis) * voigt_};
}

CauchyStress PK2Stress::to_cauchy(const Mat3& deformation_gradient) const {
    const Precision J = deformation_gradient.determinant();
    return CauchyStress::from_tensor(
        (deformation_gradient * tensor() * deformation_gradient.transpose()) / J
    );
}

PK2Stress PK2Stress::from_tensor(const Mat3& tensor) {
    return PK2Stress{Stress::from_tensor(tensor).voigt()};
}

} // namespace fem
