/**
 * @file shell_generalized_strain.h
 * @brief Declares membrane, curvature and transverse-shear strains for shells.
 *
 * The eight-component state combines three membrane strains, three curvatures
 * and two transverse shear strains. Shell elements pass this state to a
 * `ShellSection`, which returns stress resultants and the section tangent.
 *
 * @see shell_generalized_strain.cpp
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

// Section-level membrane, bending and transverse-shear strain state
struct ShellGeneralizedStrain {
    // Component order used by the internal eight-entry vector
    enum class Component : int {
        EpsilonXX = 0,
        EpsilonYY = 1,
        GammaXY   = 2,
        KappaXX   = 3,
        KappaYY   = 4,
        KappaXY   = 5,
        GammaXZ   = 6,
        GammaYZ   = 7
    };

    // Constructs a zero generalized shell strain
    ShellGeneralizedStrain() = default;

    // Constructs the state in the order defined by Component
    explicit ShellGeneralizedStrain(const Vec8& values);

    // Returns a named component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // Returns the membrane strain block [epsilon_xx, epsilon_yy, gamma_xy]
    [[nodiscard]] Vec3 membrane() const;

    // Returns the curvature block [kappa_xx, kappa_yy, kappa_xy]
    [[nodiscard]] Vec3 curvature() const;

    // Returns the transverse-shear block [gamma_xz, gamma_yz]
    [[nodiscard]] Vec2 transverse_shear() const;

    // In-plane basis transformations for membrane strain, curvature and
    // transverse-shear components. The columns of the rotation matrix contain
    // the target in-plane basis vectors expressed in the current basis.
    [[nodiscard]] ShellGeneralizedStrain transformed(const Mat2& rotation) const;
    [[nodiscard]] static Mat8            transformation(const Mat2& rotation);

    // Returns all generalized components by constant or mutable access
    [[nodiscard]] const Vec8& values() const;
    [[nodiscard]] Vec8&       values();

private:
    // [epsilon_xx, epsilon_yy, gamma_xy, kappa_xx, kappa_yy, kappa_xy, gamma_xz, gamma_yz]
    Vec8 values_{Vec8::Zero()};
};

} // namespace fem
