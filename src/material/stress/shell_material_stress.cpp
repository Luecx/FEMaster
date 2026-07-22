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

/**
 * Builds the five-component shell material stress transformation.
 *
 * The in-plane block rotates the symmetric plane-stress tensor while preserving
 * the physical-shear convention for `tau_xy`. The transverse shear entries
 * `tau_xz` and `tau_yz` are transformed as a two-dimensional vector in the
 * shell plane.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Transformation matrix for `ShellMaterialStress::values()`.
 */
Mat5 ShellMaterialStress::transformation(const Mat2& rotation) {
    const Precision a1 = rotation(0, 0);
    const Precision a2 = rotation(1, 0);
    const Precision b1 = rotation(0, 1);
    const Precision b2 = rotation(1, 1);
    const Precision two = Precision(2);

    // Rotate the in-plane physical-stress block and the transverse shear vector
    // independently
    Mat5 transformation = Mat5::Zero();
    transformation.template block<3, 3>(0, 0) <<
        a1 * a1, a2 * a2, two * a1 * a2,
        b1 * b1, b2 * b2, two * b1 * b2,
        a1 * b1, a2 * b2, a1 * b2 + a2 * b1;
    transformation.template block<2, 2>(3, 3) = rotation.transpose();

    return transformation;
}

// Expose the complete state for through-thickness integration
const Vec5& ShellMaterialStress::values() const {
    return values_;
}

Vec5& ShellMaterialStress::values() {
    return values_;
}

} // namespace fem
