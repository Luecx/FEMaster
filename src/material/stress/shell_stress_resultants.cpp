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

/**
 * Expresses the shell stress resultants in a rotated in-plane basis.
 *
 * The supplied two-by-two rotation maps current shell-basis components into the
 * target basis. Its columns contain the target in-plane basis vectors
 * expressed in the current basis.
 *
 * Membrane forces and bending moments use the physical-shear transformation for
 * symmetric in-plane stress tensors. The transverse shear forces transform as a
 * two-dimensional vector.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Stress resultants in the target basis.
 */
ShellStressResultants ShellStressResultants::transformed(const Mat2& rotation) const {
    return ShellStressResultants(transformation(rotation) * values_);
}

/**
 * Builds the eight-component shell stress-resultant transformation.
 *
 * The matrix contains two identical three-by-three physical-stress blocks for
 * membrane forces and bending moments, followed by a two-by-two vector rotation
 * for the transverse shear forces.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Transformation matrix for `ShellStressResultants::values()`.
 */
Mat8 ShellStressResultants::transformation(const Mat2& rotation) {
    const Precision a1 = rotation(0, 0);
    const Precision a2 = rotation(1, 0);
    const Precision b1 = rotation(0, 1);
    const Precision b2 = rotation(1, 1);
    const Precision two = Precision(2);

    // Build the physical-shear transformation for one symmetric in-plane
    // stress-resultant block
    Mat3 plane_stress_transform;
    plane_stress_transform <<
        a1 * a1, a2 * a2, two * a1 * a2,
        b1 * b1, b2 * b2, two * b1 * b2,
        a1 * b1, a2 * b2, a1 * b2 + a2 * b1;

    const Index membrane_start = static_cast<Index>(Component::NXX);
    const Index moment_start   = static_cast<Index>(Component::MXX);
    const Index shear_start    = static_cast<Index>(Component::QX);

    // Apply the same in-plane tensor rotation to membrane forces and moments,
    // and rotate the transverse shear force vector into the target basis
    Mat8 transformation = Mat8::Zero();
    transformation.template block<3, 3>(membrane_start, membrane_start) = plane_stress_transform;
    transformation.template block<3, 3>(moment_start, moment_start)     = plane_stress_transform;
    transformation.template block<2, 2>(shear_start, shear_start)       = rotation.transpose();

    return transformation;
}

// Expose the complete state for shell element assembly
const Vec8& ShellStressResultants::values() const {
    return values_;
}

Vec8& ShellStressResultants::values() {
    return values_;
}

} // namespace fem
