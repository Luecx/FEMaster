/**
 * @file volume_strain.h
 * @brief Declares volume Green-Lagrange strain storage and transformations.
 *
 * The material kinematics subsystem stores E = 0.5 * (F^T F - I) and its
 * increments in a six-component engineering-shear Voigt vector. This file
 * defines component access, symmetric tensor conversion, basis transformations
 * and construction from a deformation gradient. Elements and sections supply
 * the kinematics and reference bases; material models evaluate work-conjugate
 * PK2 stress and own no strain storage or element kinematics.
 *
 * @see volume_strain.cpp
 *
 * @author Finn Eggers
 * @date 09.10.2026
 */

#pragma once

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../cos/coordinate_system.h"

namespace fem {

/**
 * @brief Symmetric Green-Lagrange strain expressed in a reference basis.
 *
 * The owned Voigt vector uses [E_xx, E_yy, E_zz, 2E_yz, 2E_xz, 2E_xy].
 * Tensor conversion and orthonormal basis transformations preserve this
 * engineering-shear convention and apply equally to total strains and their
 * increments or variations. The caller distinguishes these roles; nonlinear
 * material evaluation requires a total state, while its tangent maps dE to dS.
 * Linearization about F = I recovers infinitesimal strain. Default construction
 * gives zero components without material history or external storage ownership.
 * VolumeStressPK2 supplies the work-conjugate stress representation.
 */
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

    // Computes E = 0.5 * (F^T F - I) in the reference configuration
    static VolumeStrain from_deformation_gradient(const Mat3& deformation_gradient);

    // Builds the engineering-strain Voigt transformation matrix
    static Mat6 get_transformation_matrix(const cos::Basis& from_basis,
                                          const cos::Basis& to_basis);

private:
    // [E_xx, E_yy, E_zz, 2E_yz, 2E_xz, 2E_xy], or their increments/variations
    Vec6 voigt_{Vec6::Zero()};
};

} // namespace fem
