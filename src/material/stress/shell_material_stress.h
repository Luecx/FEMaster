/**
 * @file shell_material_stress.h
 * @brief Declares the five-component stress state at a shell material point.
 *
 * The state contains two in-plane normal stresses and three shear stresses. The
 * transverse normal component is eliminated by the shell plane-stress condition.
 * Integrated shell sections evaluate this state at every point through the
 * thickness before forming force and moment resultants.
 *
 * @see shell_material_stress.cpp
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

// Local plane-stress state used at a shell material point
struct ShellMaterialStress {
    // Component order using physical shear stresses without engineering scaling
    enum class Component : int {
        XX = 0,
        YY = 1,
        XY = 2,
        XZ = 3,
        YZ = 4
    };

    // Constructs a zero shell material stress
    ShellMaterialStress() = default;

    // Constructs the state in the order defined by Component
    explicit ShellMaterialStress(const Vec5& values);

    // Returns a named component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // In-plane basis transformation for the five material-point stress
    // components. The columns of the rotation matrix contain the target
    // in-plane basis vectors expressed in the current basis.
    [[nodiscard]] static Mat5 transformation(const Mat2& rotation);

    // Returns all material-point components by constant or mutable access
    [[nodiscard]] const Vec5& values() const;
    [[nodiscard]] Vec5&       values();

protected:
    // [sigma_xx, sigma_yy, tau_xy, tau_xz, tau_yz]
    Vec5 values_{Vec5::Zero()};
};

} // namespace fem
