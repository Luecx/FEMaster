/**
 * @file shell_material_strain.h
 * @brief Declares the five-component strain state at a shell material point.
 *
 * The state contains two in-plane normal strains and the engineering shear
 * strains `gamma_xy`, `gamma_xz` and `gamma_yz`. Integrated shell sections
 * reconstruct this local state at every material point through the thickness
 * before calling a material model.
 *
 * @see shell_material_strain.cpp
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

// Local plane-stress strain state used at a shell material point
struct ShellMaterialStrain {
    // Component order, using engineering shear strains
    enum class Component : int {
        XX      = 0,
        YY      = 1,
        GammaXY = 2,
        GammaXZ = 3,
        GammaYZ = 4
    };

    // Constructs a zero shell material strain
    ShellMaterialStrain() = default;

    // Constructs the state in the order defined by Component
    explicit ShellMaterialStrain(const Vec5& values);

    // Returns a named component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // In-plane basis transformation for the five material-point strain
    // components. The columns of the rotation matrix contain the target
    // in-plane basis vectors expressed in the current basis.
    [[nodiscard]] static Mat5 transformation(const Mat2& rotation);

    // Returns all material-point components by constant or mutable access
    [[nodiscard]] const Vec5& values() const;
    [[nodiscard]] Vec5&       values();

protected:
    // [epsilon_xx, epsilon_yy, gamma_xy, gamma_xz, gamma_yz]
    Vec5 values_{Vec5::Zero()};
};

} // namespace fem
