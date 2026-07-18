/**
 * @file volume_stress_cauchy.cpp
 * @brief Implements construction and basis transformation of Cauchy stress.
 *
 * The implementation reuses `VolumeStress` storage while preserving the spatial
 * Cauchy-stress type across construction and coordinate-system transformations.
 *
 * @see volume_stress_cauchy.h
 */

#include "volume_stress_cauchy.h"

namespace fem {

// Preserve the Cauchy type for either supported representation
VolumeStressCauchy::VolumeStressCauchy(const Vec6& voigt)
    : VolumeStress(voigt) {}

VolumeStressCauchy::VolumeStressCauchy(const Mat3& tensor)
    : VolumeStress(tensor) {}

// Rotate spatial stress without losing its Cauchy-stress type
VolumeStressCauchy VolumeStressCauchy::transformed(const cos::Basis& from_basis,
                                                   const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStressCauchy(transformed);
}

} // namespace fem
