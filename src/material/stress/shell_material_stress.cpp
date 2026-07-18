/**
 * @file shell_material_stress.cpp
 * @brief Implements storage and component access for shell material stresses.
 *
 * This file maps named plane-stress components onto the five-entry vector
 * exchanged between material models and integrated shell sections.
 *
 * @see shell_material_stress.h
 */

#include "shell_material_stress.h"

namespace fem {

// Initialize the local five-component plane-stress state
ShellMaterialStress::ShellMaterialStress(const Vec5& values)
    : values_(values) {}

// Map named physical stress components onto the underlying vector
Precision& ShellMaterialStress::operator[](Component component) {
    return values_(static_cast<int>(component));
}

Precision ShellMaterialStress::operator[](Component component) const {
    return values_(static_cast<int>(component));
}

// Expose the complete state for through-thickness integration
const Vec5& ShellMaterialStress::values() const {
    return values_;
}

Vec5& ShellMaterialStress::values() {
    return values_;
}

} // namespace fem
