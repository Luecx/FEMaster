/**
 * @file shell_material_strain.h
 * @brief Declares Green-Lagrange strain components at a shell material point.
 *
 * Integrated shell sections reconstruct the five local material components
 * through the thickness from ShellGeneralizedStrain. This file provides their
 * engineering-shear storage, component access and in-plane basis transformation.
 * Material models supply the work-conjugate ShellMaterialStressPK2 response and
 * determine the missing thickness-normal strain under the plane-stress condition.
 * Generalized shell kinematics and thickness integration remain section and
 * element responsibilities.
 *
 * @see shell_material_strain.cpp
 *
 * @author Finn Eggers
 * @date 09.10.2026
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

/**
 * @brief Five Green-Lagrange strain components at a shell material point.
 *
 * Owns [E_xx, E_yy, 2E_xy, 2E_xz, 2E_yz] in the local reference material
 * basis. The missing E_zz is determined by the constitutive plane-stress
 * reduction S_zz = 0. Components may also represent increments or variations;
 * callers distinguish their role and use the material tangent for increments.
 * ShellGeneralizedStrain separately represents membrane, curvature and shear
 * quantities before reconstruction at a thickness point. Default construction
 * gives zero components, and transformation leaves the source unchanged.
 * This type owns no constitutive history and is work-conjugate to the five
 * retained components of ShellMaterialStressPK2.
 */
struct ShellMaterialStrain {
    // Component order, using engineering shear strains
    enum class Component : Index {
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

    // Expresses the Green-Lagrange components in a rotated in-plane basis
    [[nodiscard]] ShellMaterialStrain transformed(const Mat2& rotation) const;

    // Returns all material-point components by constant or mutable access
    [[nodiscard]] const Vec5& values() const;
    [[nodiscard]] Vec5&       values();

private:
    // [E_xx, E_yy, 2E_xy, 2E_xz, 2E_yz], or their increments/variations
    Vec5 values_{Vec5::Zero()};
};

} // namespace fem
