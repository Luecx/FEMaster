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

/**
 * Expresses the linearized shell material strain in a rotated in-plane basis.
 *
 * The operation preserves the linearized-strain type and applies the common
 * shell material strain transformation to the stored engineering components.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Linearized shell material strain in the target basis.
 */
ShellMaterialStrainLinearized ShellMaterialStrainLinearized::transformed(const Mat2& rotation) const {
    return ShellMaterialStrainLinearized(ShellMaterialStrain::transformation(rotation) * values_);
}

} // namespace fem
