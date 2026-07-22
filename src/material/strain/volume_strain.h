/**
 * @file volume_strain.h
 * @brief Declares three-dimensional strain storage and tensor transformations.
 *
 * `VolumeStrain` stores normal strains and engineering shear strains in a
 * six-component Voigt vector. It converts between Voigt and symmetric tensor
 * representations and transforms strains between coordinate-system bases for
 * solid material evaluation.
 *
 * @see volume_strain.cpp
 */

#pragma once

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../cos/coordinate_system.h"

namespace fem {

// Common six-component representation of a symmetric volume strain tensor
struct VolumeStrain {
    // Voigt component order, using engineering shear strains
    enum class Component : Index {
        XX      = 0,
        YY      = 1,
        ZZ      = 2,
        GammaYZ = 3,
        GammaXZ = 4,
        GammaXY = 5
    };

    // Constructs a zero volume strain
    VolumeStrain() = default;

    // Constructs from [xx, yy, zz, gamma_yz, gamma_xz, gamma_xy]
    explicit VolumeStrain(const Vec6& voigt);

    // Constructs from a symmetric tensor, converting shear entries to engineering strains
    explicit VolumeStrain(const Mat3& tensor);

    // Returns a named Voigt component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // Returns the engineering-strain Voigt vector by constant or mutable access
    [[nodiscard]] const Vec6& voigt() const;
    [[nodiscard]] Vec6&       voigt();

    // Converts the engineering-strain Voigt vector into a symmetric tensor
    [[nodiscard]] Mat3        tensor() const;

    // Expresses the strain components in another orthonormal basis
    [[nodiscard]] VolumeStrain transformed(const cos::Basis& from_basis,
                                           const cos::Basis& to_basis) const;

    // Builds the engineering-strain Voigt transformation matrix
    static Mat6 get_transformation_matrix(const cos::Basis& from_basis,
                                          const cos::Basis& to_basis);

protected:
    // [epsilon_xx, epsilon_yy, epsilon_zz, gamma_yz, gamma_xz, gamma_xy]
    Vec6 voigt_{Vec6::Zero()};
};

} // namespace fem
