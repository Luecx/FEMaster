/**
 * @file shell_material_strain_green_lagrange.cpp
 * @brief Implements construction of shell Green-Lagrange material strains.
 *
 * The implementation initializes the shared five-component shell material
 * storage while retaining the finite-strain type required for material dispatch.
 *
 * @see shell_material_strain_green_lagrange.h
 */

#include "shell_material_strain_green_lagrange.h"

namespace fem {

// Preserve the finite-strain type while using common shell material storage
ShellMaterialStrainGreenLagrange::ShellMaterialStrainGreenLagrange(const Vec5& values)
    : ShellMaterialStrain(values) {}

} // namespace fem
