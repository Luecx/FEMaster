/**
 * @file shell_material_strain.cpp
 * @brief Implements shell Green-Lagrange material strain storage and rotation.
 *
 * This material kinematics implementation provides the five-component
 * engineering-shear representation and its in-plane reference-basis rotation.
 * Total strains and their increments share this representation. Sections perform
 * thickness reconstruction and integration; constitutive models supply PK2
 * stress and enforce the missing thickness-normal plane-stress condition.
 *
 * @see shell_material_strain.h
 *
 * @author Finn Eggers
 * @date 09.10.2026
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

/**
 * @brief Expresses shell Green-Lagrange components in a rotated in-plane basis.
 *
 * The in-plane symmetric tensor and transverse shear vector are rotated using
 * the same engineering-shear convention as the stored components. The operation
 * applies equally to total strains and their increments or variations and
 * leaves the source strain and material history unchanged.
 *
 * @param rotation Target in-plane basis vectors expressed in the current basis.
 * @return Shell material strain components expressed in the target basis.
 */
ShellMaterialStrain ShellMaterialStrain::transformed(const Mat2& rotation) const {
    // Rotate the in-plane tensor components and transverse shear components.
    return ShellMaterialStrain(transformation(rotation) * values_);
}

} // namespace fem
