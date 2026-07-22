/**
 * @file shell_stress_resultants.h
 * @brief Declares membrane, bending and transverse-shear resultants for shells.
 *
 * Shell sections integrate material stresses through the thickness into three
 * membrane forces, three bending moments and two transverse shear forces. These
 * eight resultants are energetically conjugate to `ShellGeneralizedStrain`.
 *
 * @see shell_stress_resultants.cpp
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

// Generalized shell forces and moments per unit reference length
struct ShellStressResultants {
    // Component order used by the internal eight-entry vector
    enum class Component : Index {
        NXX = 0,
        NYY = 1,
        NXY = 2,
        MXX = 3,
        MYY = 4,
        MXY = 5,
        QX  = 6,
        QY  = 7
    };

    // Constructs zero shell stress resultants
    ShellStressResultants() = default;

    // Constructs resultants in the documented eight-component order
    explicit ShellStressResultants(const Vec8& values);

    // Returns a named resultant component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // Returns [Nxx, Nyy, Nxy]
    [[nodiscard]] Vec3 membrane() const;

    // Returns [Mxx, Myy, Mxy]
    [[nodiscard]] Vec3 moments() const;

    // Returns [Qx, Qy]
    [[nodiscard]] Vec2 transverse_shear() const;

    // In-plane basis transformations for membrane forces, bending moments and
    // transverse shear forces. The columns of the rotation matrix contain the
    // target in-plane basis vectors expressed in the current basis.
    [[nodiscard]] ShellStressResultants transformed(const Mat2& rotation) const;
    [[nodiscard]] static Mat8          transformation(const Mat2& rotation);

    // Returns all generalized resultants by constant or mutable access
    [[nodiscard]] const Vec8& values() const;
    [[nodiscard]] Vec8&       values();

private:
    // [Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qx, Qy]
    Vec8 values_{Vec8::Zero()};
};

} // namespace fem
