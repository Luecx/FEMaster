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

/**
 * Expresses the shell material PK2 stress in a rotated in-plane basis.
 *
 * The operation preserves the PK2-stress type and applies the common shell
 * material stress transformation to the stored physical stress components.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Shell material PK2 stress in the target basis.
 */
ShellMaterialStressPK2 ShellMaterialStressPK2::transformed(const Mat2& rotation) const {
    return ShellMaterialStressPK2(ShellMaterialStress::transformation(rotation) * values_);
}

} // namespace fem
