/**
 * @file shell_material_strain_linearized.cpp
 * @brief Implements construction of linearized shell material strains.
 *
 * The implementation initializes the shared five-component shell material
 * storage while retaining the small-strain type required for material dispatch.
 *
 * @see shell_material_strain_linearized.h
 */

#include "shell_material_strain_linearized.h"

namespace fem {

// Preserve the small-strain type while using common shell material storage
ShellMaterialStrainLinearized::ShellMaterialStrainLinearized(const Vec5& values)
    : ShellMaterialStrain(values) {}

} // namespace fem
