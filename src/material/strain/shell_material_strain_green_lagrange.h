/**
 * @file shell_material_strain_green_lagrange.h
 * @brief Declares the Green-Lagrange strain state at a shell material point.
 *
 * This type marks the five shell material components as finite-strain
 * Green-Lagrange quantities. Nonlinear integrated shell sections use it for
 * material evaluation with energetically conjugate PK2 stresses.
 *
 * @see shell_material_strain_green_lagrange.cpp
 */

#pragma once

#include "shell_material_strain.h"

namespace fem {

// Finite shell material strain used with shell PK2 stress
struct ShellMaterialStrainGreenLagrange : ShellMaterialStrain {
    // Constructs a zero Green-Lagrange shell material strain
    ShellMaterialStrainGreenLagrange() = default;

    // Constructs the finite-strain state from five engineering components
    explicit ShellMaterialStrainGreenLagrange(const Vec5& values);

    // Expresses the Green-Lagrange components in a rotated in-plane basis
    [[nodiscard]] ShellMaterialStrainGreenLagrange transformed(const Mat2& rotation) const;
};

} // namespace fem
