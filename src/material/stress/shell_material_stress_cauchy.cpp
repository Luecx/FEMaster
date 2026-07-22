/**
 * @file shell_material_stress_cauchy.cpp
 * @brief Implements construction of shell material Cauchy stress.
 *
 * The implementation initializes the shared five-component plane-stress
 * storage while retaining the spatial stress type used for material dispatch.
 *
 * @see shell_material_stress_cauchy.h
 */

#include "shell_material_stress_cauchy.h"

namespace fem {

// Preserve the Cauchy-stress type while using common shell material storage
ShellMaterialStressCauchy::ShellMaterialStressCauchy(const Vec5& values)
    : ShellMaterialStress(values) {}

/**
 * Expresses the shell material Cauchy stress in a rotated in-plane basis.
 *
 * The operation preserves the Cauchy-stress type and applies the common shell
 * material stress transformation to the stored physical stress components.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Shell material Cauchy stress in the target basis.
 */
ShellMaterialStressCauchy ShellMaterialStressCauchy::transformed(const Mat2& rotation) const {
    return ShellMaterialStressCauchy(ShellMaterialStress::transformation(rotation) * values_);
}

} // namespace fem
