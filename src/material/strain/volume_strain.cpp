/**
 * @file volume_strain.cpp
 * @brief Implements three-dimensional strain storage and tensor transformations.
 *
 * This file defines the engineering-shear Voigt convention, symmetric tensor
 * conversion and the strain transformation matrix used when solid material
 * directions differ from element coordinates.
 *
 * @see volume_strain.h
 */

#include "volume_strain.h"

namespace fem {

// Store a strain vector that already follows the engineering-shear convention
VolumeStrain::VolumeStrain(const Vec6& voigt)
    : voigt_(voigt) {}

// Convert tensor shear components into engineering shear strains
VolumeStrain::VolumeStrain(const Mat3& tensor) {
    voigt_ << tensor(0, 0),
              tensor(1, 1),
              tensor(2, 2),
              Precision(2) * tensor(1, 2),
              Precision(2) * tensor(0, 2),
              Precision(2) * tensor(0, 1);
}

// Map named components onto the underlying Voigt vector
Precision& VolumeStrain::operator[](Component component) {
    return voigt_(static_cast<int>(component));
}

Precision VolumeStrain::operator[](Component component) const {
    return voigt_(static_cast<int>(component));
}

// Expose the complete Voigt representation for material evaluation
const Vec6& VolumeStrain::voigt() const {
    return voigt_;
}

Vec6& VolumeStrain::voigt() {
    return voigt_;
}

// Recover tensor shear components by removing the engineering factor of two
Mat3 VolumeStrain::tensor() const {
    Mat3 tensor;
    tensor << voigt_(0),                  Precision(0.5) * voigt_(5), Precision(0.5) * voigt_(4),
              Precision(0.5) * voigt_(5), voigt_(1),                  Precision(0.5) * voigt_(3),
              Precision(0.5) * voigt_(4), Precision(0.5) * voigt_(3), voigt_(2);
    return tensor;
}

// Rotate the strain representation between orthonormal bases
VolumeStrain VolumeStrain::transformed(const cos::Basis& from_basis,
                                       const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStrain(transformed);
}

// Build the Voigt transformation for engineering shear strains
Mat6 VolumeStrain::get_transformation_matrix(const cos::Basis& from_basis,
                                             const cos::Basis& to_basis) {
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

    Mat6 transformation;
    transformation <<
        R11 * R11, R12 * R12, R13 * R13, R12 * R13, R11 * R13, R11 * R12,
        R21 * R21, R22 * R22, R23 * R23, R22 * R23, R21 * R23, R21 * R22,
        R31 * R31, R32 * R32, R33 * R33, R32 * R33, R31 * R33, R31 * R32,
        Precision(2) * R21 * R31, Precision(2) * R22 * R32, Precision(2) * R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
        Precision(2) * R11 * R31, Precision(2) * R12 * R32, Precision(2) * R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
        Precision(2) * R11 * R21, Precision(2) * R12 * R22, Precision(2) * R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;
    return transformation;
}

} // namespace fem
