/**
 * @file shell_generalized_strain.cpp
 * @brief Implements storage and component access for generalized shell strains.
 *
 * The implementation maps named membrane, curvature and transverse-shear
 * components onto the section-level eight-entry vector used by shell assembly.
 *
 * @see shell_generalized_strain.h
 */

#include "shell_generalized_strain.h"

namespace fem {

// Initialize membrane, curvature and transverse-shear components
ShellGeneralizedStrain::ShellGeneralizedStrain(const Vec8& values)
    : values_(values) {}

// Map named components onto the underlying vector
Precision& ShellGeneralizedStrain::operator[](Component component) {
    return values_(static_cast<int>(component));
}

Precision ShellGeneralizedStrain::operator[](Component component) const {
    return values_(static_cast<int>(component));
}

// Extract the three section-kinematic blocks
Vec3 ShellGeneralizedStrain::membrane() const {
    return values_.template segment<3>(0);
}

Vec3 ShellGeneralizedStrain::curvature() const {
    return values_.template segment<3>(3);
}

Vec2 ShellGeneralizedStrain::transverse_shear() const {
    return values_.template segment<2>(6);
}

// Expose the complete state for ShellSection evaluation
const Vec8& ShellGeneralizedStrain::values() const {
    return values_;
}

Vec8& ShellGeneralizedStrain::values() {
    return values_;
}

} // namespace fem
