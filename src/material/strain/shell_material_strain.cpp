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

/**
 * Builds the five-component shell material strain transformation.
 *
 * The in-plane block rotates the symmetric plane-stress strain tensor while
 * preserving the engineering-shear convention for `gamma_xy`. The transverse
 * shear entries `gamma_xz` and `gamma_yz` are transformed as a two-dimensional
 * vector in the shell plane.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Transformation matrix for `ShellMaterialStrain::values()`.
 */
Mat5 ShellMaterialStrain::transformation(const Mat2& rotation) {
    const Precision a1 = rotation(0, 0);
    const Precision a2 = rotation(1, 0);
    const Precision b1 = rotation(0, 1);
    const Precision b2 = rotation(1, 1);
    const Precision two = Precision(2);

    // Rotate the in-plane engineering-strain block and the transverse shear
    // vector independently
    Mat5 transformation = Mat5::Zero();
    transformation.template block<3, 3>(0, 0) <<
        a1 * a1,       a2 * a2,       a1 * a2,
        b1 * b1,       b2 * b2,       b1 * b2,
        two * a1 * b1, two * a2 * b2, a1 * b2 + a2 * b1;
    transformation.template block<2, 2>(3, 3) = rotation.transpose();

    return transformation;
}

// Expose the complete state for material evaluation
const Vec5& ShellMaterialStrain::values() const {
    return values_;
}

Vec5& ShellMaterialStrain::values() {
    return values_;
}

} // namespace fem
