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

// DEBUG A/B: force the constitutive Green-Lagrange strain to zero while the
// element still builds B_GL from the actual deformation gradient.
VolumeStrainGreenLagrange VolumeStrainGreenLagrange::from_deformation_gradient(
    const Mat3& deformation_gradient
) {
    (void)deformation_gradient;
    return VolumeStrainGreenLagrange{};
}

} // namespace fem
