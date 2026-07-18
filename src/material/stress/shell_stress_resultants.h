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
    // Constructs zero shell stress resultants
    ShellStressResultants() = default;

    // Constructs resultants in the documented eight-component order
    explicit ShellStressResultants(const Vec8& values);

    // Returns [Nxx, Nyy, Nxy]
    [[nodiscard]] Vec3 membrane() const;

    // Returns [Mxx, Myy, Mxy]
    [[nodiscard]] Vec3 moments() const;

    // Returns [Qx, Qy]
    [[nodiscard]] Vec2 transverse_shear() const;

    // Returns all generalized resultants by constant or mutable access
    [[nodiscard]] const Vec8& values() const;
    [[nodiscard]] Vec8&       values();

private:
    // [Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qx, Qy]
    Vec8 values_{Vec8::Zero()};
};

} // namespace fem
