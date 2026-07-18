/**
 * @file volume_stress_cauchy.h
 * @brief Declares the three-dimensional Cauchy stress measure.
 *
 * Cauchy stress is defined in the current configuration and represents true
 * force per current area. Linearized solid materials return it directly, while
 * finite-strain output can obtain it by pushing forward PK2 stress.
 *
 * @see volume_stress_cauchy.cpp
 */

#pragma once

#include "volume_stress.h"

namespace fem {

// Spatial true stress defined in the current configuration
struct VolumeStressCauchy : VolumeStress {
    // Constructs a zero Cauchy stress
    VolumeStressCauchy() = default;

    // Constructs from the physical-shear Voigt representation
    explicit VolumeStressCauchy(const Vec6& voigt);

    // Constructs from the symmetric Cauchy stress tensor
    explicit VolumeStressCauchy(const Mat3& tensor);

    // Expresses the Cauchy stress in another orthonormal basis
    [[nodiscard]] VolumeStressCauchy transformed(const cos::Basis& from_basis,
                                                 const cos::Basis& to_basis) const;
};

} // namespace fem
