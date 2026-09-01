/**
 * @file section_solid.cpp
 * @brief Implements oriented constitutive evaluation for solid sections.
 *
 * Global linearized or Green-Lagrange strains are transformed into the optional
 * reference material basis, evaluated through the assigned elasticity model and
 * returned as global Cauchy or second Piola-Kirchhoff stress with a consistently
 * transformed tangent. The implementation also differentiates the transformed
 * linear tangent with respect to additional material-rotation parameters.
 *
 * Material-point history remains owned by the active nonlinear solution. The
 * section receives separate selected input and output rows and forwards both to
 * the constitutive model.
 *
 * @see SolidSection
 * @see material::Elasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "section_solid.h"

#include "../core/logging.h"

#include <sstream>

namespace fem {

/**
 * Returns the reference material basis at one physical point.
 *
 * Spatial coordinate systems are evaluated at the point expressed in their own
 * local coordinates. Without an assigned orientation, the global Cartesian
 * basis is the material basis.
 *
 * @param position_reference Physical reference position of the material point.
 * @return Material basis whose columns are expressed in global coordinates.
 */
Mat3 SolidSection::section_orientation_basis(const Vec3& position_reference) const {
    if (!orientation_) {
        return Mat3::Identity();
    }

    const Vec3 point_local = orientation_->to_local(position_reference);
    return orientation_->get_axes(point_local);
}

/**
 * Evaluates linearized solid stress and tangent in global coordinates.
 *
 * The optional section orientation and additional element rotation define the
 * current material axes in the reference configuration. Engineering strain is
 * transformed from global to material coordinates, evaluated as Cauchy stress
 * through the assigned elasticity and transformed back to global coordinates.
 * The tangent follows the identical chain of stress and strain transformations.
 *
 * The selected old material-point state row is passed read-only to the
 * constitutive law, which writes updated history to the corresponding new row.
 *
 * @param position_reference Physical reference position of the material point.
 * @param additional_rotation Element-provided rotation applied after the section basis.
 * @param strain_global Linearized engineering strain in global coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Material-point output state row.
 * @param stress_global Cauchy stress returned in global coordinates.
 * @param tangent_global Consistent global tangent mapping global strain to stress.
 */
void SolidSection::evaluate(const Vec3&                   position_reference,
                            const Mat3&                   additional_rotation,
                            const VolumeStrainLinearized& strain_global,
                            const Precision*              old_state,
                            Precision*                    new_state,
                            VolumeStressCauchy&           stress_global,
                            Mat6&                         tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");

    auto elasticity = material_->elasticity();
    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const Vec6                   strain_material_values = strain_transform * strain_global.voigt();
    const VolumeStrainLinearized strain_material(strain_material_values);

    VolumeStressCauchy stress_material;
    Mat6               tangent_material;
    elasticity->evaluate(strain_material, old_state, new_state, stress_material, tangent_material);

    const Vec6 stress_global_values = stress_transform * stress_material.voigt();
    stress_global  = VolumeStressCauchy(stress_global_values);
    tangent_global = stress_transform * tangent_material * strain_transform;
}

/**
 * Evaluates Total-Lagrangian solid stress and tangent in global coordinates.
 *
 * DEBUG A/B: retain the old GL element kinematics, but evaluate the constitutive
 * tangent through the linearized/Cauchy material overload. This isolates whether
 * the shell_solid_tie regression comes from the material strain/stress API rather
 * than the Green-Lagrange B matrix or current/reference geometry.
 */
void SolidSection::evaluate(const Vec3&                      position_reference,
                            const Mat3&                      additional_rotation,
                            const VolumeStrainGreenLagrange& strain_global,
                            const Precision*                 old_state,
                            Precision*                       new_state,
                            VolumeStressPK2&                 stress_global,
                            Mat6&                            tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");

    auto elasticity = material_->elasticity();
    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const Vec6 strain_material_values = strain_transform * strain_global.voigt();
    const VolumeStrainLinearized strain_material(strain_material_values);

    VolumeStressCauchy stress_material;
    Mat6               tangent_material;

    (void)new_state;
    elasticity->evaluate(strain_material, old_state, nullptr, stress_material, tangent_material);

    const Vec6 stress_global_values = stress_transform * stress_material.voigt();
    stress_global  = VolumeStressPK2(stress_global_values);
    tangent_global = stress_transform * tangent_material * strain_transform;
}

/**
 * As stated in `SolidElement<N>::compute_compliance_angle_derivative`, computes
 * the derivative of the material tangent with respect to the three additional
 * material-orientation angles.
 *
 * @param position_reference Reference position of the material point.
 * @param additional_rotation Current additional material rotation.
 * @param additional_rotation_derivatives Derivatives of the additional rotation.
 * @param old_state Immutable material-point input state selected by the element.
 * @param new_state Material-point output state selected by the element.
 * @return Derivatives of the global tangent with respect to the three angles.
 */
std::array<Mat6, 3> SolidSection::tangent_rotation_derivatives(
    const Vec3&                position_reference,
    const Mat3&                additional_rotation,
    const std::array<Mat3, 3>& additional_rotation_derivatives,
    const Precision*           old_state,
    Precision*                 new_state
) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");

    auto elasticity = material_->elasticity();

    VolumeStrainLinearized zero_strain;
    VolumeStressCauchy     zero_stress;
    Mat6                   tangent_material;
    elasticity->evaluate(zero_strain, old_state, new_state, zero_stress, tangent_material);

    const Mat3 section_basis  = section_orientation_basis(position_reference);
    const Mat3 material_basis = section_basis * additional_rotation;

    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    auto strain_transform_derivative = [&](const Mat3& material_basis_derivative) {
        Mat6 derivative;

        for (Index component = 0; component < 6; ++component) {
            Vec6 unit = Vec6::Zero();
            unit(component) = Precision(1);

            const Mat3 global_tensor = VolumeStrain(unit).tensor();
            const Mat3 local_derivative =
                material_basis_derivative.transpose() * global_tensor * material_basis
                + material_basis.transpose() * global_tensor * material_basis_derivative;

            derivative.col(component) = VolumeStrain(local_derivative).voigt();
        }

        return derivative;
    };

    auto stress_transform_derivative = [&](const Mat3& material_basis_derivative) {
        Mat6 derivative;

        for (Index component = 0; component < 6; ++component) {
            Vec6 unit       = Vec6::Zero();
            unit(component) = Precision(1);

            const Mat3 local_tensor = VolumeStress(unit).tensor();
            const Mat3 global_derivative
                = material_basis_derivative * local_tensor * material_basis.transpose()
                + material_basis            * local_tensor * material_basis_derivative.transpose();

            derivative.col(component) = VolumeStress(global_derivative).voigt();
        }

        return derivative;
    };

    std::array<Mat6, 3> tangent_derivatives;

    for (Index i = 0; i < 3; ++i) {
        const Mat3 material_basis_derivative = section_basis * additional_rotation_derivatives[i];
        const Mat6 strain_derivative         = strain_transform_derivative(material_basis_derivative);
        const Mat6 stress_derivative         = stress_transform_derivative(material_basis_derivative);

        tangent_derivatives[i] =
            stress_derivative * tangent_material * strain_transform
            + stress_transform * tangent_material * strain_derivative;
    }

    return tangent_derivatives;
}

void SolidSection::info() {
    logging::info(true, "SolidSection:");
    logging::info(true, "   Material   : ", (material_    ? material_   ->name : "-"));
    logging::info(true, "   Region     : ", (region_      ? region_     ->name : "-"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
}

std::string SolidSection::str() const {
    std::ostringstream os;

    os << "SolidSection: material=" << (material_    ? material_   ->name : std::string("-"))
       << ", orientation="          << (orientation_ ? orientation_->name : std::string("-"))
       << ", region="               << (region_      ? region_->name : std::string("-"))
       << " ("                      << (region_      ? static_cast<int>(region_->size()) : 0) << ")";

    return os.str();
}
} // namespace fem
