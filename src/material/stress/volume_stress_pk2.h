/**
 * @file volume_stress_pk2.h
 * @brief Declares the three-dimensional second Piola-Kirchhoff stress measure.
 *
 * PK2 stress is defined in the reference configuration and is energetically
 * conjugate to Green-Lagrange strain. The type can be transformed between
 * material bases or pushed forward to Cauchy stress for spatial output.
 *
 * @see volume_stress_pk2.cpp
 */

#pragma once

#include "volume_stress.h"
#include "volume_stress_cauchy.h"

namespace fem {

// Material stress defined in the reference configuration
struct VolumeStressPK2 : VolumeStress {
    // Constructs a zero PK2 stress
    VolumeStressPK2() = default;

    // Constructs from the physical-shear Voigt representation
    explicit VolumeStressPK2(const Vec6& voigt);

    // Constructs from the symmetric PK2 stress tensor
    explicit VolumeStressPK2(const Mat3& tensor);

    // Expresses the PK2 stress in another reference basis
    [[nodiscard]] VolumeStressPK2 transformed(const cos::Basis& from_basis,
                                              const cos::Basis& to_basis) const;

    // Pushes PK2 stress forward using sigma = F * S * F^T / det(F)
    [[nodiscard]] VolumeStressCauchy to_cauchy(const Mat3& deformation_gradient) const;
};

} // namespace fem
