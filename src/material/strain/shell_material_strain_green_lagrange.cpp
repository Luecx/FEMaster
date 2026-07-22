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

/**
 * Expresses the Green-Lagrange shell material strain in a rotated in-plane
 * basis.
 *
 * The operation preserves the Green-Lagrange type and applies the common shell
 * material strain transformation to the stored engineering components.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Green-Lagrange shell material strain in the target basis.
 */
ShellMaterialStrainGreenLagrange ShellMaterialStrainGreenLagrange::transformed(const Mat2& rotation) const {
    return ShellMaterialStrainGreenLagrange(ShellMaterialStrain::transformation(rotation) * values_);
}

} // namespace fem
