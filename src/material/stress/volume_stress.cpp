/**
 * @file volume_stress.cpp
 * @brief Implements three-dimensional stress storage and tensor transformations.
 *
 * This file defines the physical-shear Voigt convention, symmetric tensor
 * conversion and stress transformation matrix used by solid materials and
 * stress-output paths.
 *
 * @see volume_stress.h
 */

#include "volume_stress.h"

namespace fem {

// Store a stress vector that already follows the physical-shear convention
VolumeStress::VolumeStress(const Vec6& voigt)
    : voigt_(voigt) {}

// Extract the six independent components of a symmetric stress tensor
VolumeStress::VolumeStress(const Mat3& tensor) {
    voigt_ << tensor(0, 0),
              tensor(1, 1),
              tensor(2, 2),
              tensor(1, 2),
              tensor(0, 2),
              tensor(0, 1);
}

// Map named components onto the underlying Voigt vector
Precision& VolumeStress::operator[](Component component) {
    return voigt_(static_cast<int>(component));
}

Precision VolumeStress::operator[](Component component) const {
    return voigt_(static_cast<int>(component));
}

// Expose the complete Voigt representation for material and output paths
const Vec6& VolumeStress::voigt() const {
    return voigt_;
}

Vec6& VolumeStress::voigt() {
    return voigt_;
}

// Reconstruct the symmetric tensor without engineering-shear scaling
Mat3 VolumeStress::tensor() const {
    Mat3 tensor;
    tensor << voigt_(0), voigt_(5), voigt_(4),
              voigt_(5), voigt_(1), voigt_(3),
              voigt_(4), voigt_(3), voigt_(2);
    return tensor;
}

// Rotate the stress representation between orthonormal bases
VolumeStress VolumeStress::transformed(const cos::Basis& from_basis,
                                       const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStress(transformed);
}

// Build the Voigt transformation for physical shear stresses
Mat6 VolumeStress::get_transformation_matrix(const cos::Basis& from_basis,
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
        R11 * R11, R12 * R12, R13 * R13, Precision(2) * R12 * R13, Precision(2) * R11 * R13, Precision(2) * R11 * R12,
        R21 * R21, R22 * R22, R23 * R23, Precision(2) * R22 * R23, Precision(2) * R21 * R23, Precision(2) * R21 * R22,
        R31 * R31, R32 * R32, R33 * R33, Precision(2) * R32 * R33, Precision(2) * R31 * R33, Precision(2) * R31 * R32,
        R21 * R31, R22 * R32, R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
        R11 * R31, R12 * R32, R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
        R11 * R21, R12 * R22, R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;
    return transformation;
}

} // namespace fem
