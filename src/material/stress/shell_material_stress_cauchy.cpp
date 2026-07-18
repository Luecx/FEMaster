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

} // namespace fem
