/**
 * @file volume_strain_green_lagrange.cpp
 * @brief Implements Green-Lagrange strain construction and transformations.
 *
 * The implementation constructs the strain tensor from a deformation gradient
 * and preserves the Green-Lagrange type when transforming between bases.
 *
 * @see volume_strain_green_lagrange.h
 */

#include "volume_strain_green_lagrange.h"

namespace fem {

// Preserve the Green-Lagrange type for either supported representation
VolumeStrainGreenLagrange::VolumeStrainGreenLagrange(const Vec6& voigt)
    : VolumeStrain(voigt) {}

VolumeStrainGreenLagrange::VolumeStrainGreenLagrange(const Mat3& tensor)
    : VolumeStrain(tensor) {}

// Rotate the material strain without losing its finite-strain type
VolumeStrainGreenLagrange VolumeStrainGreenLagrange::transformed(const cos::Basis& from_basis,
                                                                 const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStrainGreenLagrange(transformed);
}

// Form E = 0.5 * (F^T F - I) in the reference configuration
VolumeStrainGreenLagrange VolumeStrainGreenLagrange::from_deformation_gradient(
    const Mat3& deformation_gradient
) {
    const Mat3 strain_tensor = Precision(0.5)
        * (deformation_gradient.transpose() * deformation_gradient - Mat3::Identity());
    return VolumeStrainGreenLagrange(strain_tensor);
}

} // namespace fem
