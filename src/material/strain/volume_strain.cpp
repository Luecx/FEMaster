/**
 * @file volume_strain.cpp
 * @brief Implements volume Green-Lagrange strain storage and transformations.
 *
 * This material kinematics implementation defines engineering-shear Voigt/tensor
 * conversion, reference-basis transformations and E = 0.5 * (F^T F - I).
 * The same linear conversions apply to strain increments and variations.
 * Elements and sections supply deformation gradients and bases; constitutive
 * evaluation and material history are outside this type's responsibility.
 *
 * @see volume_strain.h
 *
 * @author Finn Eggers
 * @date 09.10.2026
 */

#include "volume_strain.h"

namespace fem {

// Store a strain vector that already follows the engineering-shear convention
VolumeStrain::VolumeStrain(const Vec6& voigt)
    : voigt_(voigt) {}

/**
 * @brief Converts a symmetric Green-Lagrange tensor to engineering Voigt form.
 *
 * Stores the three normal entries and twice the yz, xz and xy tensor entries.
 * The supplied tensor represents a total strain or an increment/variation in
 * the caller's reference basis; symmetry is a caller precondition.
 *
 * @param tensor Symmetric Green-Lagrange tensor or its increment/variation.
 */
VolumeStrain::VolumeStrain(const Mat3& tensor) {
    // Double shear entries so the strain vector pairs with physical PK2 shear stress.
    voigt_ << tensor(0, 0),
              tensor(1, 1),
              tensor(2, 2),
              Precision(2) * tensor(1, 2),
              Precision(2) * tensor(0, 2),
              Precision(2) * tensor(0, 1);
}

// Map named components onto the underlying Voigt vector
Precision& VolumeStrain::operator[](Component component) {
    return voigt_(static_cast<int>(component));
}

Precision VolumeStrain::operator[](Component component) const {
    return voigt_(static_cast<int>(component));
}

// Expose the complete Voigt representation for material evaluation
const Vec6& VolumeStrain::voigt() const {
    return voigt_;
}

Vec6& VolumeStrain::voigt() {
    return voigt_;
}

/**
 * @brief Recovers the symmetric strain tensor in the stored reference basis.
 *
 * Normal components are unchanged; engineering shear components are halved.
 * Applies to total Green-Lagrange strains and their increments or variations.
 *
 * @return Symmetric tensor without modifying the stored components.
 */
Mat3 VolumeStrain::tensor() const {
    // Undo the engineering-shear factor of two in both symmetric tensor entries.
    Mat3 tensor;
    tensor << voigt_(0),                  Precision(0.5) * voigt_(5), Precision(0.5) * voigt_(4),
              Precision(0.5) * voigt_(5), voigt_(1),                  Precision(0.5) * voigt_(3),
              Precision(0.5) * voigt_(4), Precision(0.5) * voigt_(3), voigt_(2);
    return tensor;
}

/**
 * @brief Expresses Green-Lagrange components in another reference basis.
 *
 * For R = Q_to^T Q_from, the tensor transforms as E_to = R E_from R^T.
 * The Voigt operator preserves engineering shear and also transforms increments.
 *
 * @param from_basis Orthonormal source basis expressed in global coordinates.
 * @param to_basis Orthonormal target basis expressed in global coordinates.
 * @return Transformed strain, leaving the source components unchanged.
 */
VolumeStrain VolumeStrain::transformed(const cos::Basis& from_basis,
                                       const cos::Basis& to_basis) const {
    // Apply the reference-basis rotation in engineering Voigt representation.
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStrain(transformed);
}

/**
 * @brief Builds the engineering-strain operator for reference-basis rotation.
 *
 * Expands E_to = R E_from R^T with R = Q_to^T Q_from into six-component
 * Voigt form. Factors of two in the shear rows account for stored engineering
 * components; no constitutive state or strain components are modified.
 *
 * @param from_basis Orthonormal source basis expressed in global coordinates.
 * @param to_basis Orthonormal target basis expressed in global coordinates.
 * @return Operator mapping source to target engineering-strain components.
 */
Mat6 VolumeStrain::get_transformation_matrix(const cos::Basis& from_basis,
                                             const cos::Basis& to_basis) {
    // Resolve source basis vectors in the target reference basis.
    const Mat3 R = to_basis.transpose() * from_basis;

    const Precision R11 = R(0, 0);
    const Precision R12 = R(0, 1);
    const Precision R13 = R(0, 2);
    const Precision R21 = R(1, 0);
    const Precision R22 = R(1, 1);
    const Precision R23 = R(1, 2);
    const Precision R31 = R(2, 0);
    const Precision R32 = R(2, 1);
    const Precision R33 = R(2, 2);

    // Expand the symmetric tensor rotation, including engineering-shear factors.
    Mat6 transformation;
    transformation <<
        R11 * R11, R12 * R12, R13 * R13, R12 * R13, R11 * R13, R11 * R12,
        R21 * R21, R22 * R22, R23 * R23, R22 * R23, R21 * R23, R21 * R22,
        R31 * R31, R32 * R32, R33 * R33, R32 * R33, R31 * R33, R31 * R32,
        Precision(2) * R21 * R31, Precision(2) * R22 * R32, Precision(2) * R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
        Precision(2) * R11 * R31, Precision(2) * R12 * R32, Precision(2) * R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
        Precision(2) * R11 * R21, Precision(2) * R12 * R22, Precision(2) * R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;
    return transformation;
}

/**
 * @brief Constructs Green-Lagrange strain from a deformation gradient.
 *
 * F maps reference line elements to the current configuration. The right
 * Cauchy-Green tensor C = F^T F therefore gives E = 0.5 * (C - I) in the
 * reference basis. The constructor converts its shear entries to engineering
 * components. Element kinematics remain responsible for validating F.
 *
 * @param deformation_gradient Deformation gradient in the supplied bases.
 * @return Total Green-Lagrange strain expressed in the reference basis.
 */
VolumeStrain VolumeStrain::from_deformation_gradient(const Mat3& deformation_gradient) {
    // Form the reference metric and subtract the undeformed metric.
    const Mat3 strain_tensor = Precision(0.5)
        * (deformation_gradient.transpose() * deformation_gradient - Mat3::Identity());
    return VolumeStrain(strain_tensor);
}

} // namespace fem
