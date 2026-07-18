/**
 * @file volume_stress_pk2.cpp
 * @brief Implements PK2 stress transformations and conversion to Cauchy stress.
 *
 * PK2 components can be rotated between reference bases and pushed forward with
 * a deformation gradient for output in the current configuration.
 *
 * @see volume_stress_pk2.h
 */

#include "volume_stress_pk2.h"

#include <Eigen/LU>

namespace fem {

// Preserve the PK2 type for either supported representation
VolumeStressPK2::VolumeStressPK2(const Vec6& voigt)
    : VolumeStress(voigt) {}

VolumeStressPK2::VolumeStressPK2(const Mat3& tensor)
    : VolumeStress(tensor) {}

// Rotate material stress without losing its PK2-stress type
VolumeStressPK2 VolumeStressPK2::transformed(const cos::Basis& from_basis,
                                             const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStressPK2(transformed);
}

// Push reference stress forward with sigma = F * S * F^T / det(F)
VolumeStressCauchy VolumeStressPK2::to_cauchy(const Mat3& deformation_gradient) const {
    const Precision determinant = deformation_gradient.determinant();
    const Mat3 cauchy =
        deformation_gradient * tensor() * deformation_gradient.transpose() / determinant;
    return VolumeStressCauchy(cauchy);
}

} // namespace fem
