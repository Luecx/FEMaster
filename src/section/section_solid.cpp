/**
 * @file section_solid.cpp
 * @brief Implements oriented constitutive updates and recovery for solid sections.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#include "section_solid.h"

#include "../core/logging.h"

#include <sstream>

namespace fem {

Mat3 SolidSection::section_orientation_basis(const Vec3& position_reference) const {
    if (!orientation_) return Mat3::Identity();
    const Vec3 point_local = orientation_->to_local(position_reference);
    return orientation_->get_axes(point_local);
}

void SolidSection::update(const Vec3& position_reference,
                          const Mat3& additional_rotation,
                          const VolumeStrainLinearized& strain_global,
                          const Precision* state_old,
                          Precision* state_new,
                          VolumeStressCauchy& stress_global,
                          Mat6& tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    const auto& law = material_->constitutive_law();
    logging::error(law.supports_volume_linearized(),
        "SolidSection constitutive law does not support linearized volume evaluation");

    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());
    const VolumeStrainLinearized strain_material(strain_transform * strain_global.voigt());

    VolumeStressCauchy stress_material;
    Mat6 tangent_material;
    law.update(strain_material, state_old, state_new, stress_material, tangent_material);

    stress_global = VolumeStressCauchy(stress_transform * stress_material.voigt());
    tangent_global = stress_transform * tangent_material * strain_transform;
}

void SolidSection::update(const Vec3& position_reference,
                          const Mat3& additional_rotation,
                          const VolumeStrainGreenLagrange& strain_global,
                          const Precision* state_old,
                          Precision* state_new,
                          VolumeStressPK2& stress_global,
                          Mat6& tangent_global) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    const auto& law = material_->constitutive_law();
    logging::error(law.supports_volume_green_lagrange(),
        "SolidSection constitutive law does not support Green-Lagrange volume evaluation");

    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());
    const VolumeStrainGreenLagrange strain_material(strain_transform * strain_global.voigt());

    VolumeStressPK2 stress_material;
    Mat6 tangent_material;
    law.update(strain_material, state_old, state_new, stress_material, tangent_material);

    stress_global = VolumeStressPK2(stress_transform * stress_material.voigt());
    tangent_global = stress_transform * tangent_material * strain_transform;
}

void SolidSection::recover(const Vec3& position_reference,
                           const Mat3& additional_rotation,
                           const VolumeStrainLinearized& strain_global,
                           const Precision* state,
                           VolumeStressCauchy& stress_global) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    const auto& law = material_->constitutive_law();
    logging::error(law.supports_volume_linearized(),
        "SolidSection constitutive law does not support linearized volume recovery");

    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());
    const VolumeStrainLinearized strain_material(strain_transform * strain_global.voigt());

    VolumeStressCauchy stress_material;
    law.recover(strain_material, state, stress_material);
    stress_global = VolumeStressCauchy(stress_transform * stress_material.voigt());
}

void SolidSection::recover(const Vec3& position_reference,
                           const Mat3& additional_rotation,
                           const VolumeStrainGreenLagrange& strain_global,
                           const Precision* state,
                           VolumeStressPK2& stress_global) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    const auto& law = material_->constitutive_law();
    logging::error(law.supports_volume_green_lagrange(),
        "SolidSection constitutive law does not support Green-Lagrange volume recovery");

    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());
    const VolumeStrainGreenLagrange strain_material(strain_transform * strain_global.voigt());

    VolumeStressPK2 stress_material;
    law.recover(strain_material, state, stress_material);
    stress_global = VolumeStressPK2(stress_transform * stress_material.voigt());
}

Mat6 SolidSection::elastic_tangent_reference(const Vec3& position_reference,
                                             const Mat3& additional_rotation) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    auto elasticity = material_->elastic_backbone();
    logging::error(elasticity->supports_volume_linearized(),
        "SolidSection elastic backbone does not support linearized volume evaluation");
    logging::error(elasticity->state_size() == 0,
        "SolidSection elastic reference tangent requires a stateless elastic backbone");

    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    VolumeStrainLinearized zero_strain;
    VolumeStressCauchy zero_stress;
    Mat6 tangent_material;
    elasticity->evaluate(zero_strain, nullptr, zero_stress, tangent_material);
    return stress_transform * tangent_material * strain_transform;
}

std::array<Mat6, 3> SolidSection::tangent_rotation_derivatives(
    const Vec3& position_reference,
    const Mat3& additional_rotation,
    const std::array<Mat3, 3>& additional_rotation_derivatives
) const {
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(!material_->has_plasticity(),
        "Material-orientation tangent derivatives are currently elastic-only");
    auto elasticity = material_->elastic_backbone();
    logging::error(elasticity->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");
    logging::error(elasticity->state_size() == 0,
        "Material-orientation tangent derivatives require stateless elasticity");

    VolumeStrainLinearized zero_strain;
    VolumeStressCauchy zero_stress;
    Mat6 tangent_material;
    elasticity->evaluate(zero_strain, nullptr, zero_stress, tangent_material);

    const Mat3 section_basis = section_orientation_basis(position_reference);
    const Mat3 material_basis = section_basis * additional_rotation;
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    auto strain_transform_derivative = [&](const Mat3& dQ) {
        Mat6 derivative;
        for (Index component = 0; component < 6; ++component) {
            Vec6 unit = Vec6::Zero();
            unit(component) = Precision(1);
            const Mat3 global_tensor = VolumeStrain(unit).tensor();
            const Mat3 local_derivative = dQ.transpose() * global_tensor * material_basis
                                        + material_basis.transpose() * global_tensor * dQ;
            derivative.col(component) = VolumeStrain(local_derivative).voigt();
        }
        return derivative;
    };

    auto stress_transform_derivative = [&](const Mat3& dQ) {
        Mat6 derivative;
        for (Index component = 0; component < 6; ++component) {
            Vec6 unit = Vec6::Zero();
            unit(component) = Precision(1);
            const Mat3 local_tensor = VolumeStress(unit).tensor();
            const Mat3 global_derivative = dQ * local_tensor * material_basis.transpose()
                                         + material_basis * local_tensor * dQ.transpose();
            derivative.col(component) = VolumeStress(global_derivative).voigt();
        }
        return derivative;
    };

    std::array<Mat6, 3> derivatives;
    for (Index i = 0; i < 3; ++i) {
        const Mat3 dQ = section_basis * additional_rotation_derivatives[i];
        derivatives[i] = stress_transform_derivative(dQ) * tangent_material * strain_transform
                       + stress_transform * tangent_material * strain_transform_derivative(dQ);
    }
    return derivatives;
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
