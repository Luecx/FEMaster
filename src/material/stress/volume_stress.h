/**
 * @file volume_stress.h
 * @brief Declares three-dimensional stress storage and tensor transformations.
 *
 * `VolumeStress` stores a symmetric stress tensor in six-component Voigt form.
 * Unlike strain, shear stresses are physical tensor components and are not
 * multiplied by two. Solid material and output paths use the tensor conversion
 * and basis transformation facilities provided here.
 *
 * @see volume_stress.cpp
 */

#pragma once

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../cos/coordinate_system.h"

namespace fem {

// Common six-component representation of a symmetric volume stress tensor
struct VolumeStress {
    // Voigt component order using physical shear stresses
    enum class Component : int {
        XX = 0,
        YY = 1,
        ZZ = 2,
        YZ = 3,
        XZ = 4,
        XY = 5
    };

    // Constructs a zero volume stress
    VolumeStress() = default;

    // Constructs from [xx, yy, zz, yz, xz, xy]
    explicit VolumeStress(const Vec6& voigt);

    // Constructs from a symmetric stress tensor
    explicit VolumeStress(const Mat3& tensor);

    // Returns a named Voigt component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // Returns the stress Voigt vector by constant or mutable access
    [[nodiscard]] const Vec6& voigt() const;
    [[nodiscard]] Vec6&       voigt();

    // Converts the stress Voigt vector into a symmetric tensor
    [[nodiscard]] Mat3        tensor() const;

    // Expresses the stress components in another orthonormal basis
    [[nodiscard]] VolumeStress transformed(const cos::Basis& from_basis,
                                            const cos::Basis& to_basis) const;

    // Builds the physical-shear stress transformation matrix
    static Mat6 get_transformation_matrix(const cos::Basis& from_basis,
                                          const cos::Basis& to_basis);

protected:
    // [sigma_xx, sigma_yy, sigma_zz, tau_yz, tau_xz, tau_xy]
    Vec6 voigt_{Vec6::Zero()};
};

} // namespace fem
