/**
 * @file shell_stress_resultants.cpp
 * @brief Implements storage and block access for generalized shell resultants.
 *
 * The implementation exposes membrane, moment and transverse-shear blocks from
 * the eight-entry state assembled by shell sections.
 *
 * @see shell_stress_resultants.h
 */

#include "shell_stress_resultants.h"

namespace fem {

// Initialize all generalized shell forces and moments
ShellStressResultants::ShellStressResultants(const Vec8& values)
    : values_(values) {}

// Map named components onto the underlying vector
Precision& ShellStressResultants::operator[](Component component) {
    return values_(static_cast<int>(component));
}

Precision ShellStressResultants::operator[](Component component) const {
    return values_(static_cast<int>(component));
}

// Extract the three physical resultant blocks
Vec3 ShellStressResultants::membrane() const {
    return values_.template segment<3>(static_cast<int>(Component::NXX));
}

Vec3 ShellStressResultants::moments() const {
    return values_.template segment<3>(static_cast<int>(Component::MXX));
}

Vec2 ShellStressResultants::transverse_shear() const {
    return values_.template segment<2>(static_cast<int>(Component::QX));
}

// Expose the complete state for shell element assembly
const Vec8& ShellStressResultants::values() const {
    return values_;
}

Vec8& ShellStressResultants::values() {
    return values_;
}

} // namespace fem
