/**
 * @file section_solid.cpp
 * @brief Implements oriented constitutive evaluation for solid sections.
 *
 * Read-only evaluation and history integration share the same material-basis
 * transformations. Only `integrate()` may write constitutive history, and it
 * forwards immutable source and separate target rows to the material law.
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

Mat3 SolidSection::section_orientation_basis(const Vec3& position_reference) const {
    if (!orientation_) {
        return Mat3::Identity();
    }

    const Vec3 point_local = orientation_->to_local(position_reference);
    return orientation_->get_axes(point_local);
}

void SolidSection::evaluate(const Vec3&                   position_reference,
                            const Mat3&                   additional_rotation,
                            const VolumeStrainLinearized& strain_global,
                            const Precision*              state,
                            VolumeStressCauchy&           stress_global,
                            Mat6&                         tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
                   "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
                   "SolidSection material does not support linearized volume evaluation");

    const auto elasticity = material_->elasticity();
    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const VolumeStrainLinearized strain_material(strain_transform * strain_global.voigt());
    VolumeStressCauchy stress_material;
    Mat6 tangent_material;

    elasticity->evaluate(strain_material, state, stress_material, tangent_material);

    stress_global = VolumeStressCauchy(stress_transform * stress_material.voigt());
    tangent_global = stress_transform * tangent_material * strain_transform;
}

void SolidSection::evaluate(const Vec3&                      position_reference,
                            const Mat3&                      additional_rotation,
                            const VolumeStrainGreenLagrange& strain_global,
                            const Precision*                 state,
                            VolumeStressPK2&                 stress_global,
                            Mat6&                            tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
                   "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_green_lagrange(),
                   "SolidSection material does not support Green-Lagrange volume evaluation");

    const auto elasticity = material_->elasticity();
    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const VolumeStrainGreenLagrange strain_material(strain_transform * strain_global.voigt());
    VolumeStressPK2 stress_material;
    Mat6 tangent_material;

    elasticity->evaluate(strain_material, state, stress_material, tangent_material);

    stress_global = VolumeStressPK2(stress_transform * stress_material.voigt());
    tangent_global = stress_transform * tangent_material * strain_transform;
}

void SolidSection::integrate(const Vec3&                   position_reference,
                             const Mat3&                   additional_rotation,
                             const VolumeStrainLinearized& strain_global,
                             const Precision*              state,
                             Precision*                    target_state,
                             VolumeStressCauchy&           stress_global,
                             Mat6&                         tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
                   "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
                   "SolidSection material does not support linearized volume integration");

    const auto elasticity = material_->elasticity();
    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const VolumeStrainLinearized strain_material(strain_transform * strain_global.voigt());
    VolumeStressCauchy stress_material;
    Mat6 tangent_material;

    elasticity->integrate(
        strain_material, state, target_state, stress_material, tangent_material
    );

    stress_global = VolumeStressCauchy(stress_transform * stress_material.voigt());
    tangent_global = stress_transform * tangent_material * strain_transform;
}

void SolidSection::integrate(const Vec3&                      position_reference,
                             const Mat3&                      additional_rotation,
                             const VolumeStrainGreenLagrange& strain_global,
                             const Precision*                 state,
                             Precision*                       target_state,
                             VolumeStressPK2&                 stress_global,
                             Mat6&                            tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
                   "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_green_lagrange(),
                   "SolidSection material does not support Green-Lagrange volume integration");

    const auto elasticity = material_->elasticity();
    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const VolumeStrainGreenLagrange strain_material(strain_transform * strain_global.voigt());
    VolumeStressPK2 stress_material;
    Mat6 tangent_material;

    elasticity->integrate(
        strain_material, state, target_state, stress_material, tangent_material
    );

    stress_global = VolumeStressPK2(stress_transform * stress_material.voigt());
    tangent_global = stress_transform * tangent_material * strain_transform;
}

std::array<Mat6, 3> SolidSection::tangent_rotation_derivatives(
    const Vec3&                position_reference,
    const Mat3&                additional_rotation,
    const std::array<Mat3, 3>& additional_rotation_derivatives,
    const Precision*           state
) const {
    logging::error(material_ && material_->has_elasticity(),
                   "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
                   "SolidSection material does not support linearized volume evaluation");

    const auto elasticity = material_->elasticity();

    VolumeStrainLinearized zero_strain;
    VolumeStressCauchy zero_stress;
    Mat6 tangent_material;
    elasticity->evaluate(zero_strain, state, zero_stress, tangent_material);

    const Mat3 section_basis  = section_orientation_basis(position_reference);
    const Mat3 material_basis = section_basis * additional_rotation;
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

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
            Vec6 unit = Vec6::Zero();
            unit(component) = Precision(1);

            const Mat3 local_tensor = VolumeStress(unit).tensor();
            const Mat3 global_derivative =
                material_basis_derivative * local_tensor * material_basis.transpose()
              + material_basis * local_tensor * material_basis_derivative.transpose();

            derivative.col(component) = VolumeStress(global_derivative).voigt();
        }
        return derivative;
    };

    std::array<Mat6, 3> tangent_derivatives;
    for (Index i = 0; i < 3; ++i) {
        const Mat3 material_basis_derivative = section_basis * additional_rotation_derivatives[i];
        const Mat6 strain_derivative = strain_transform_derivative(material_basis_derivative);
        const Mat6 stress_derivative = stress_transform_derivative(material_basis_derivative);

        tangent_derivatives[i] =
            stress_derivative * tangent_material * strain_transform
          + stress_transform * tangent_material * strain_derivative;
    }

    return tangent_derivatives;
}

void SolidSection::info() {
    logging::info(true, "SolidSection:");
    logging::info(true, "   Material   : ", (material_ ? material_->name : "-"));
    logging::info(true, "   Region     : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
}

std::string SolidSection::str() const {
    std::ostringstream os;
    os << "SolidSection: material=" << (material_ ? material_->name : std::string("-"))
       << ", orientation=" << (orientation_ ? orientation_->name : std::string("-"))
       << ", region=" << (region_ ? region_->name : std::string("-"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")";
    return os.str();
}

} // namespace fem
