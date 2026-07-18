/**
 * @file shell_material_stress_pk2.cpp
 * @brief Implements construction of shell material PK2 stress.
 *
 * The implementation initializes the shared five-component plane-stress
 * storage while retaining the reference-configuration stress type.
 *
 * @see shell_material_stress_pk2.h
 */

#include "shell_material_stress_pk2.h"

namespace fem {

// Preserve the PK2-stress type while using common shell material storage
ShellMaterialStressPK2::ShellMaterialStressPK2(const Vec5& values)
    : ShellMaterialStress(values) {}

} // namespace fem
