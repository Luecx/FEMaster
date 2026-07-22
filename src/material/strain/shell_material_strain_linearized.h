/**
 * @file shell_material_strain_linearized.h
 * @brief Declares the linearized strain state at a shell material point.
 *
 * This type marks the local shell components as small-strain quantities. Linear
 * integrated shell sections use it when requesting Cauchy stress and a
 * linearized material tangent.
 *
 * @see shell_material_strain_linearized.cpp
 */

#pragma once

#include "shell_material_strain.h"

namespace fem {

// Small shell material strain used with shell Cauchy stress
struct ShellMaterialStrainLinearized : ShellMaterialStrain {
    // Constructs a zero linearized shell material strain
    ShellMaterialStrainLinearized() = default;

    // Constructs the small-strain state from five engineering components
    explicit ShellMaterialStrainLinearized(const Vec5& values);

    // Expresses the small-strain components in a rotated in-plane basis
    [[nodiscard]] ShellMaterialStrainLinearized transformed(const Mat2& rotation) const;
};

} // namespace fem
