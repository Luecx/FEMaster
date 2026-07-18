/**
 * @file shell_material_strain.cpp
 * @brief Implements storage and component access for shell material strains.
 *
 * This file maps named shell material components onto the five-entry vector
 * exchanged between integrated shell sections and material models.
 *
 * @see shell_material_strain.h
 */

#include "shell_material_strain.h"

namespace fem {

// Initialize the local five-component plane-stress state
ShellMaterialStrain::ShellMaterialStrain(const Vec5& values)
    : values_(values) {}

// Map named components onto the underlying vector
Precision& ShellMaterialStrain::operator[](Component component) {
    return values_(static_cast<int>(component));
}

Precision ShellMaterialStrain::operator[](Component component) const {
    return values_(static_cast<int>(component));
}

// Expose the complete state for material evaluation
const Vec5& ShellMaterialStrain::values() const {
    return values_;
}

Vec5& ShellMaterialStrain::values() {
    return values_;
}

} // namespace fem
