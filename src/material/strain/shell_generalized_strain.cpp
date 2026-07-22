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
    return values_.template segment<3>(static_cast<int>(Component::EpsilonXX));
}

Vec3 ShellGeneralizedStrain::curvature() const {
    return values_.template segment<3>(static_cast<int>(Component::KappaXX));
}

Vec2 ShellGeneralizedStrain::transverse_shear() const {
    return values_.template segment<2>(static_cast<int>(Component::GammaXZ));
}

/**
 * Expresses the generalized shell strain in a rotated in-plane basis.
 *
 * The supplied two-by-two rotation maps current shell-basis components into the
 * target basis. Its columns contain the target in-plane basis vectors
 * expressed in the current basis.
 *
 * Membrane strain and curvature use the engineering-shear transformation for
 * symmetric in-plane tensors. The transverse shear block is transformed as a
 * two-dimensional vector.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Generalized strain components in the target basis.
 */
ShellGeneralizedStrain ShellGeneralizedStrain::transformed(const Mat2& rotation) const {
    return ShellGeneralizedStrain(transformation(rotation) * values_);
}

/**
 * Builds the eight-component generalized shell strain transformation.
 *
 * The matrix contains two identical three-by-three engineering-strain blocks
 * for membrane strain and curvature, followed by a two-by-two vector rotation
 * for the transverse engineering shear strains.
 *
 * @param rotation In-plane target basis vectors expressed in the current basis.
 * @return Transformation matrix for `ShellGeneralizedStrain::values()`.
 */
Mat8 ShellGeneralizedStrain::transformation(const Mat2& rotation) {
    const Precision a1  = rotation(0, 0);
    const Precision a2  = rotation(1, 0);
    const Precision b1  = rotation(0, 1);
    const Precision b2  = rotation(1, 1);
    const Precision two = Precision(2);

    // Build the engineering-shear transformation for one symmetric in-plane
    // strain block
    Mat3 plane_strain_transform;
    plane_strain_transform <<
        a1 * a1,       a2 * a2,       a1 * a2,
        b1 * b1,       b2 * b2,       b1 * b2,
        two * a1 * b1, two * a2 * b2, a1 * b2 + a2 * b1;

    const Index membrane_start  = static_cast<Index>(Component::EpsilonXX);
    const Index curvature_start = static_cast<Index>(Component::KappaXX);
    const Index shear_start     = static_cast<Index>(Component::GammaXZ);

    // Apply the same in-plane tensor rotation to membrane strain and curvature,
    // and rotate the transverse shear vector into the target basis
    Mat8 transformation = Mat8::Zero();
    transformation.template block<3, 3>(membrane_start, membrane_start)   = plane_strain_transform;
    transformation.template block<3, 3>(curvature_start, curvature_start) = plane_strain_transform;
    transformation.template block<2, 2>(shear_start, shear_start)         = rotation.transpose();

    return transformation;
}

// Expose the complete state for ShellSection evaluation
const Vec8& ShellGeneralizedStrain::values() const {
    return values_;
}

Vec8& ShellGeneralizedStrain::values() {
    return values_;
}

} // namespace fem
