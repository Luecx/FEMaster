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

// Extract the three physical resultant blocks
Vec3 ShellStressResultants::membrane() const {
    return values_.template segment<3>(0);
}

Vec3 ShellStressResultants::moments() const {
    return values_.template segment<3>(3);
}

Vec2 ShellStressResultants::transverse_shear() const {
    return values_.template segment<2>(6);
}

// Expose the complete state for shell element assembly
const Vec8& ShellStressResultants::values() const {
    return values_;
}

Vec8& ShellStressResultants::values() {
    return values_;
}

} // namespace fem
