/**
 * @file volume_strain_linearized.cpp
 * @brief Implements construction and transformations of linearized volume strains.
 *
 * The implementation initializes the shared Voigt storage and preserves the
 * linearized strain type when material orientations require a basis change.
 *
 * @see volume_strain_linearized.h
 */

#include "volume_strain_linearized.h"

namespace fem {

// Preserve the linearized type for either supported representation
VolumeStrainLinearized::VolumeStrainLinearized(const Vec6& voigt)
    : VolumeStrain(voigt) {}

VolumeStrainLinearized::VolumeStrainLinearized(const Mat3& tensor)
    : VolumeStrain(tensor) {}

// Rotate the small-strain state without losing its material-dispatch type
VolumeStrainLinearized VolumeStrainLinearized::transformed(const cos::Basis& from_basis,
                                                           const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStrainLinearized(transformed);
}

} // namespace fem
