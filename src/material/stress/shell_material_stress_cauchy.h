/**
 * @file shell_material_stress_cauchy.h
 * @brief Declares the Cauchy stress state at a shell material point.
 *
 * This type identifies the five plane-stress components as spatial true
 * stresses. Linearized integrated shell sections use it with the corresponding
 * linearized shell material strain.
 *
 * @see shell_material_stress_cauchy.cpp
 */

#pragma once

#include "shell_material_stress.h"

namespace fem {

// Shell material Cauchy stress defined in the current configuration
struct ShellMaterialStressCauchy : ShellMaterialStress {
    // Constructs a zero shell material Cauchy stress
    ShellMaterialStressCauchy() = default;

    // Constructs the Cauchy stress from five plane-stress components
    explicit ShellMaterialStressCauchy(const Vec5& values);

    // Expresses the Cauchy stress components in a rotated in-plane basis
    [[nodiscard]] ShellMaterialStressCauchy transformed(const Mat2& rotation) const;
};

} // namespace fem
