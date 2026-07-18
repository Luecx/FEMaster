/**
 * @file volume_strain_linearized.h
 * @brief Declares the three-dimensional linearized strain measure.
 *
 * This type identifies the common six-component volume strain as a small-strain
 * quantity. Linear solid elements use it for material evaluation with Cauchy
 * stress and a linearized tangent.
 *
 * @see volume_strain_linearized.cpp
 */

#pragma once

#include "volume_strain.h"

namespace fem {

// Linearized volume strain used with Cauchy stress
struct VolumeStrainLinearized : VolumeStrain {
    // Constructs a zero linearized volume strain
    VolumeStrainLinearized() = default;

    // Constructs from the engineering-strain Voigt representation
    explicit VolumeStrainLinearized(const Vec6& voigt);

    // Constructs from the symmetric linearized strain tensor
    explicit VolumeStrainLinearized(const Mat3& tensor);

    // Expresses the linearized strain in another orthonormal basis
    [[nodiscard]] VolumeStrainLinearized transformed(const cos::Basis& from_basis,
                                                      const cos::Basis& to_basis) const;
};

} // namespace fem
