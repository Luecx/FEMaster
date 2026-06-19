/**
 * @file strain.cpp
 * @brief Implements strain rotation utilities in Voigt notation.
 *
 * @see src/material/strain.h
 * @see src/cos/coordinate_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "strain.h"

namespace fem {

Strain::Strain(const Vec6& voigt)
    : voigt_(voigt) {}

Strain::Strain(const Mat3& tensor) {
    voigt_ << tensor(0, 0),
              tensor(1, 1),
              tensor(2, 2),
              Precision(2) * tensor(1, 2),
              Precision(2) * tensor(2, 0),
              Precision(2) * tensor(0, 1);
}

Precision& Strain::operator[](Component component) {
    return voigt_(static_cast<int>(component));
}

Precision Strain::operator[](Component component) const {
    return voigt_(static_cast<int>(component));
}

Precision Strain::operator[](ComponentDerived component) const {
    switch (component) {
        case ComponentDerived::TrueYZ:
            return Precision(0.5) * voigt_(3);
        case ComponentDerived::TrueZX:
            return Precision(0.5) * voigt_(4);
        case ComponentDerived::TrueXY:
            return Precision(0.5) * voigt_(5);
    }
    return Precision(0);
}

void Strain::set(ComponentDerived component, Precision value) {
    switch (component) {
        case ComponentDerived::TrueYZ:
            voigt_(3) = Precision(2) * value;
            return;
        case ComponentDerived::TrueZX:
            voigt_(4) = Precision(2) * value;
            return;
        case ComponentDerived::TrueXY:
            voigt_(5) = Precision(2) * value;
            return;
    }
}

const Vec6& Strain::voigt() const {
    return voigt_;
}

Vec6& Strain::voigt() {
    return voigt_;
}

Mat3 Strain::tensor() const {
    Mat3 tensor;
    tensor << voigt_(0),                  Precision(0.5) * voigt_(5), Precision(0.5) * voigt_(4),
              Precision(0.5) * voigt_(5), voigt_(1),                  Precision(0.5) * voigt_(3),
              Precision(0.5) * voigt_(4), Precision(0.5) * voigt_(3), voigt_(2);
    return tensor;
}

Strain Strain::transformed(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) const {
    const Vec6 transformed = Strain::get_transformation_matrix(from_basis, to_basis) * voigt_;
    return Strain(transformed);
}

Mat6 Strain::get_transformation_matrix(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) {
    const Mat3 R = to_basis.transpose() * from_basis;

    Precision R11 = R(0, 0), R12 = R(0, 1), R13 = R(0, 2);
    Precision R21 = R(1, 0), R22 = R(1, 1), R23 = R(1, 2);
    Precision R31 = R(2, 0), R32 = R(2, 1), R33 = R(2, 2);

    StaticMatrix<6, 6> T_eps;
    T_eps <<
        R11 * R11, R12 * R12, R13 * R13, R12 * R13, R11 * R13, R11 * R12,
        R21 * R21, R22 * R22, R23 * R23, R22 * R23, R21 * R23, R21 * R22,
        R31 * R31, R32 * R32, R33 * R33, R32 * R33, R31 * R33, R31 * R32,
        2 * R21 * R31, 2 * R22 * R32, 2 * R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
        2 * R11 * R31, 2 * R12 * R32, 2 * R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
        2 * R11 * R21, 2 * R12 * R22, 2 * R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;

    return T_eps;
}

Mat6 Strain::get_transformation_matrix(const cos::Basis& basis) {
    Precision R11 = basis(0, 0), R12 = basis(0, 1), R13 = basis(0, 2);
    Precision R21 = basis(1, 0), R22 = basis(1, 1), R23 = basis(1, 2);
    Precision R31 = basis(2, 0), R32 = basis(2, 1), R33 = basis(2, 2);

    StaticMatrix<6, 6> T_eps;
    T_eps <<
        R11 * R11, R12 * R12, R13 * R13, R12 * R13, R11 * R13, R11 * R12,
        R21 * R21, R22 * R22, R23 * R23, R22 * R23, R21 * R23, R21 * R22,
        R31 * R31, R32 * R32, R33 * R33, R32 * R33, R31 * R33, R31 * R32,
        2 * R21 * R31, 2 * R22 * R32, 2 * R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
        2 * R11 * R31, 2 * R12 * R32, 2 * R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
        2 * R11 * R21, 2 * R12 * R22, 2 * R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;

    return T_eps;
}

LinearizedStrain::LinearizedStrain(const Vec6& voigt)
    : Strain(voigt) {}

LinearizedStrain::LinearizedStrain(const Mat3& tensor)
    : Strain(tensor) {}

LinearizedStrain LinearizedStrain::transformed(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) const {
    const Vec6 transformed = Strain::get_transformation_matrix(from_basis, to_basis) * voigt_;
    return LinearizedStrain(transformed);
}

GreenLagrangeStrain::GreenLagrangeStrain(const Vec6& voigt)
    : Strain(voigt) {}

GreenLagrangeStrain::GreenLagrangeStrain(const Mat3& tensor)
    : Strain(tensor) {}

GreenLagrangeStrain GreenLagrangeStrain::transformed(
    const cos::Basis& from_basis,
    const cos::Basis& to_basis
) const {
    const Vec6 transformed = Strain::get_transformation_matrix(from_basis, to_basis) * voigt_;
    return GreenLagrangeStrain(transformed);
}

GreenLagrangeStrain GreenLagrangeStrain::from_deformation_gradient(
    const Mat3& deformation_gradient
) {
    const Mat3 strain_tensor = Precision(0.5)
        * (deformation_gradient.transpose() * deformation_gradient - Mat3::Identity());
    return GreenLagrangeStrain(strain_tensor);
}
} // namespace fem
