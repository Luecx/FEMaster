/**
 * @file section_solid.cpp
 * @brief Implements solid section reporting.
 *
 * @see src/section/section_solid.h
 * @author Finn Eggers
 * @date 28.04.2026
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
                            VolumeStressCauchy&            stress_global,
                            Mat6&                          tangent_global) const {

    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");
    logging::error(material_->elasticity()->state_size() == 0           ,
        "SolidSection does not yet provide integration-point material state storage");

    // get the elastic material model from the material
    auto elasticity = material_->elasticity();

    // since the given strains are in global coordinates but we may have a local material orientation, we need to
    // get the material basis. we add in the additional rotation
    const Mat3 material_basis   = section_orientation_basis(position_reference) * additional_rotation;

    // technically we can use strain_global.transformed() but we need the transformation later again
    // so strain_transform can transform the strain to the material system and stress_transform can undo that.
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    // actually transform strain to local material basis
    const Vec6                   strain_material_values = strain_transform * strain_global.voigt();
    const VolumeStrainLinearized strain_material(strain_material_values);

    // prepare cauchy stress and tangent as output for the elasticity
    VolumeStressCauchy stress_material;
    Mat6               tangent_material;
    // call into material model
    elasticity->evaluate(
        strain_material,
        nullptr,
        nullptr,
        stress_material,
        tangent_material
    );

    // turn the stress back into global coordinates
    const Vec6 stress_global_values = stress_transform * stress_material.voigt();
    stress_global  = VolumeStressCauchy(stress_global_values);
    // Transform the material tangent from the local material basis to the global basis:
    // C_global = T_stress * C_material * T_strain
    tangent_global = stress_transform * tangent_material * strain_transform;
}
void SolidSection::evaluate(const Vec3&                      position_reference,
                            const Mat3&                      additional_rotation,
                            const VolumeStrainGreenLagrange& strain_global,
                            VolumeStressPK2&                 stress_global,
                            Mat6&                            tangent_global) const {

    logging::error(material_ && material_->has_elasticity(),
                   "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_green_lagrange(),
        "SolidSection material does not support Green-Lagrange volume evaluation");
    logging::error(material_->elasticity()->state_size() == 0,
        "SolidSection does not yet provide integration-point material state storage");

    // get the elastic material model from the material
    auto elasticity = material_->elasticity();

    // since the given strains are in global coordinates but we may have a local material orientation, we need to
    // get the material basis. we add in the additional rotation
    const Mat3 material_basis = section_orientation_basis(position_reference) * additional_rotation;

    // technically we can use strain_global.transformed() but we need the transformation later again
    // so strain_transform can transform the strain to the material system and stress_transform can undo that.
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    // actually transform Green-Lagrange strain to local material basis
    const Vec6                      strain_material_values = strain_transform * strain_global.voigt();
    const VolumeStrainGreenLagrange strain_material(strain_material_values);

    // prepare PK2 stress and tangent as output for the elasticity
    VolumeStressPK2 stress_material;
    Mat6            tangent_material;
    // call into material model
    elasticity->evaluate(
        strain_material,
        nullptr,
        nullptr,
        stress_material,
        tangent_material
    );

    // turn the stress back into global coordinates
    const Vec6 stress_global_values = stress_transform * stress_material.voigt();
    stress_global  = VolumeStressPK2(stress_global_values);
    // Transform the material tangent from the local material basis to the global basis:
    // C_global = T_stress * C_material * T_strain
    tangent_global = stress_transform * tangent_material * strain_transform;
}

/**
 *
 * As stated in @SolidElement<N>::compute_compliance_angle_derivative, we need the derivative of the tangent stiffness
 * w.r.t to the 3 angles.
 *
 *
 * @param position_reference
 * @param additional_rotation
 * @param additional_rotation_derivatives
 * @return
 */
std::array<Mat6, 3> SolidSection::tangent_rotation_derivatives(
    const Vec3&                position_reference,
    const Mat3&                additional_rotation,
    const std::array<Mat3, 3>& additional_rotation_derivatives
) const {
    // a few checks
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");
    logging::error(material_->elasticity()->state_size() == 0,
        "SolidSection does not yet provide integration-point material state storage");

    // get the elastic material model from the material
    auto elasticity = material_->elasticity();

    // The global material tangent is obtained by transforming the local tangent:
    //
    //     C_global = T_stress(Q) * C_material * T_strain(Q)
    //
    // Here, Q is the material basis expressed in global coordinates. This function assumes that
    // C_material is constant with respect to the rotation parameters, as is the case for linear
    // elasticity. Therefore, it is sufficient to evaluate the material tangent at zero strain.
    VolumeStrainLinearized zero_strain;
    VolumeStressCauchy     zero_stress;
    Mat6                   tangent_material;
    elasticity->evaluate(
        zero_strain,
        nullptr,
        nullptr,
        zero_stress,
        tangent_material
    );

    // The material basis is composed of the prescribed section orientation Q_section and an
    // additional rotation R:
    //
    //     Q = Q_section * R
    //
    // Since Q_section is independent of the additional rotation parameters q_i, its derivative is:
    //
    //     dQ/dq_i = Q_section * dR/dq_i
    const Mat3 section_basis  = section_orientation_basis(position_reference);
    const Mat3 material_basis = section_basis * additional_rotation;

    // T_strain transforms a global strain into the local material basis:
    //
    //     eps_material = Q^T * eps_global * Q
    //                  = T_strain(Q) * eps_global
    //
    // T_stress transforms a local material stress back into the global basis:
    //
    //     sigma_global = Q * sigma_material * Q^T
    //                  = T_stress(Q) * sigma_material
    const Mat6 strain_transform = VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform = VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    // Compute the derivative dT_strain/dq for a given material-basis derivative dQ/dq.
    //
    // Starting from:
    //
    //     eps_material = Q^T * eps_global * Q
    //
    // and applying the product rule gives:
    //
    //     d(eps_material)/dq
    //         = (dQ/dq)^T * eps_global * Q
    //         + Q^T * eps_global * dQ/dq
    //
    // Since T_strain is a 6x6 linear operator, its derivative is constructed column by column.
    // Each Voigt unit vector represents one basis strain, and the transformed derivative of this
    // basis strain forms the corresponding column of dT_strain/dq.
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

    // Compute the derivative dT_stress/dq for a given material-basis derivative dQ/dq.
    //
    // Starting from:
    //
    //     sigma_global = Q * sigma_material * Q^T
    //
    // and applying the product rule gives:
    //
    //     d(sigma_global)/dq
    //         = dQ/dq * sigma_material * Q^T
    //         + Q * sigma_material * (dQ/dq)^T
    //
    // As for the strain transformation, the derivative of the 6x6 stress transformation is
    // constructed column by column using the six Voigt basis stresses.
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

    // Differentiate the globally transformed material tangent:
    //
    //     C_global = T_stress * C_material * T_strain
    //
    // Assuming dC_material/dq_i = 0, the product rule gives:
    //
    //     dC_global/dq_i
    //         = dT_stress/dq_i * C_material * T_strain
    //         + T_stress * C_material * dT_strain/dq_i
    //
    // The three returned matrices therefore contain the derivatives of C_global with respect to
    // the three supplied additional-rotation parameters.
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
       << ", region="               << (region_      ? region_     ->name : std::string("-"))
       << " ("                      << (region_      ? static_cast<int>(region_->size()) : 0) << ")";

    return os.str();
}
} // namespace fem
