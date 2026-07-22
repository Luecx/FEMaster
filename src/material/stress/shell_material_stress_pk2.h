/**
 * @file shell_material_stress_pk2.h
 * @brief Declares the PK2 stress state at a shell material point.
 *
 * This type identifies the five plane-stress components as reference-configuration
 * stresses energetically conjugate to shell Green-Lagrange material strains.
 * Nonlinear integrated shell sections use it for finite-strain materials.
 *
 * @see shell_material_stress_pk2.cpp
 */

#pragma once

#include "shell_material_stress.h"

namespace fem {

// Shell material PK2 stress defined in the reference configuration
struct ShellMaterialStressPK2 : ShellMaterialStress {
    // Constructs a zero shell material PK2 stress
    ShellMaterialStressPK2() = default;

    // Constructs the PK2 stress from five plane-stress components
    explicit ShellMaterialStressPK2(const Vec5& values);

    // Expresses the PK2 stress components in a rotated in-plane basis
    [[nodiscard]] ShellMaterialStressPK2 transformed(const Mat2& rotation) const;
};

} // namespace fem
