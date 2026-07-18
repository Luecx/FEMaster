/**
 * @file volume_strain_green_lagrange.h
 * @brief Declares the three-dimensional Green-Lagrange strain measure.
 *
 * The type represents `E = 0.5 * (F^T F - I)` in the common six-component
 * volume-strain storage. Geometrically nonlinear solid elements use it for
 * finite-strain material evaluation with PK2 stress.
 *
 * @see volume_strain_green_lagrange.cpp
 */

#pragma once

#include "volume_strain.h"

namespace fem {

// Green-Lagrange volume strain energetically conjugate to PK2 stress
struct VolumeStrainGreenLagrange : VolumeStrain {
    // Constructs a zero Green-Lagrange strain
    VolumeStrainGreenLagrange() = default;

    // Constructs from the engineering-strain Voigt representation
    explicit VolumeStrainGreenLagrange(const Vec6& voigt);

    // Constructs from the symmetric Green-Lagrange strain tensor
    explicit VolumeStrainGreenLagrange(const Mat3& tensor);

    // Expresses the Green-Lagrange strain in another orthonormal basis
    [[nodiscard]] VolumeStrainGreenLagrange transformed(const cos::Basis& from_basis,
                                                        const cos::Basis& to_basis) const;

    // Computes E = 0.5 * (F^T F - I) from a deformation gradient
    static VolumeStrainGreenLagrange from_deformation_gradient(const Mat3& deformation_gradient);
};

} // namespace fem
