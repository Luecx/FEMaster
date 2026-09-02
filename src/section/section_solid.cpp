/**
 * @file section_solid.cpp
 * @brief Implements oriented constitutive evaluation for solid sections.
 *
 * Global linearized or Green-Lagrange strains are transformed into the optional
 * reference material basis, evaluated through the assigned elasticity model and
 * returned as global Cauchy or second Piola-Kirchhoff stress. Consistent tangents
 * follow the same transformation path only when requested by the caller.
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
 * Evaluates linearized solid response in global coordinates.
 *
 * The optional section orientation and additional element rotation define the
 * material basis in the reference configuration. Engineering strain is
 * transformed from global to material coordinates and evaluated as Cauchy stress.
 * The material stress is then transformed back into global coordinates.
 *
 * If a tangent is requested, the material derivative is transformed by
 *
 *     C_global = T_stress C_material T_strain.
 *
 * A null tangent pointer is forwarded to the elasticity model, so a stress-only
 * query does not construct or transform a material tangent.
 *
 * @param position_reference Physical reference position of the material point.
 * @param additional_rotation Element-provided rotation applied after the section basis.
 * @param strain_global Linearized engineering strain in global coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Optional material-point output state row.
 * @param stress_global Cauchy stress returned in global coordinates.
 * @param tangent_global Optional consistent global tangent.
 */
void SolidSection::evaluate(const Vec3&                   position_reference,
                            const Mat3&                   additional_rotation,
                            const VolumeStrainLinearized& strain_global,
                            const Precision*              old_state,
                            Precision*                    new_state,
                            VolumeStressCauchy&           stress_global,
                            Mat6*                         tangent_global) const {
    // Validate the requested constitutive formulation before performing any
    // coordinate transformation.
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");

    auto elasticity = material_->elasticity();

    // Compose the user-defined section orientation with the element-provided
    // additional material rotation.
    const Mat3 material_basis =
        section_orientation_basis(position_reference) * additional_rotation;

    // Transform global engineering strain into the material basis. The stress
    // transformation is always required because stress is mandatory output.
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const VolumeStrainLinearized strain_material(
        strain_transform * strain_global.voigt()
    );

    // Evaluate the constitutive response. Material tangent storage exists only
    // when the caller actually requested a global tangent.
    VolumeStressCauchy stress_material;
    Mat6               tangent_material;

    elasticity->evaluate(
        strain_material,
        old_state,
        new_state,
        stress_material,
        tangent_global != nullptr ? &tangent_material : nullptr
    );

    // Transform the physical Cauchy stress back into global coordinates.
    stress_global = VolumeStressCauchy(
        stress_transform * stress_material.voigt()
    );

    if (tangent_global != nullptr) {
        // Apply the complete material-to-global chain rule only when required:
        //
        //     d sigma_g / d epsilon_g
        //         = T_stress C_material T_strain.
        *tangent_global = stress_transform * tangent_material * strain_transform;
    }
}

/**
 * Evaluates Total-Lagrangian solid response in global reference coordinates.
 *
 * Green-Lagrange strain is transformed into the material basis, evaluated by the
 * assigned elasticity and returned as second Piola-Kirchhoff stress. Both strain
 * and stress remain reference-configuration quantities, so the constitutive
 * tangent is transformed by the same engineering-Voigt chain rule used in the
 * linearized case.
 *
 * A null tangent pointer is propagated to the material model and skips the final
 * tangent transformation. The stress and constitutive state update remain
 * identical to a full tangent evaluation.
 *
 * @param position_reference Physical reference position of the material point.
 * @param additional_rotation Element-provided rotation applied after the section basis.
 * @param strain_global Green-Lagrange strain in global reference coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Optional material-point output state row.
 * @param stress_global PK2 stress returned in global reference coordinates.
 * @param tangent_global Optional consistent global material tangent `dS/dE`.
 */
void SolidSection::evaluate(const Vec3&                      position_reference,
                            const Mat3&                      additional_rotation,
                            const VolumeStrainGreenLagrange& strain_global,
                            const Precision*                 old_state,
                            Precision*                       new_state,
                            VolumeStressPK2&                 stress_global,
                            Mat6*                            tangent_global) const {
    // Validate material availability and the exact finite-strain formulation.
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_green_lagrange(),
        "SolidSection material does not support Green-Lagrange volume evaluation");

    auto elasticity = material_->elasticity();

    // Build the reference material basis and both engineering-Voigt transforms.
    const Mat3 material_basis =
        section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    const VolumeStrainGreenLagrange strain_material(
        strain_transform * strain_global.voigt()
    );

    // Forward tangent optionality all the way into the constitutive model.
    VolumeStressPK2 stress_material;
    Mat6            tangent_material;

    elasticity->evaluate(
        strain_material,
        old_state,
        new_state,
        stress_material,
        tangent_global != nullptr ? &tangent_material : nullptr
    );

    // PK2 stress is transformed between reference material bases without changing
    // its stress measure.
    stress_global = VolumeStressPK2(
        stress_transform * stress_material.voigt()
    );

    if (tangent_global != nullptr) {
        *tangent_global = stress_transform * tangent_material * strain_transform;
    }
}

/**
 * Differentiates the globally transformed linear material tangent with respect
 * to three supplied additional material-orientation parameters.
 *
 * The global tangent is
 *
 *     C_global = T_stress(Q) C_material T_strain(Q).
 *
 * The linear material operator is independent of the orientation parameters, so
 * differentiation gives
 *
 *     dC_global/dq
 *         = dT_stress/dq C_material T_strain
 *         + T_stress C_material dT_strain/dq.
 *
 * This sensitivity explicitly requires `C_material`; unlike ordinary stress-only
 * constitutive calls, a tangent is therefore requested from the elasticity at
 * zero strain.
 *
 * @param position_reference Reference position of the material point.
 * @param additional_rotation Current additional material rotation.
 * @param additional_rotation_derivatives Derivatives of the additional rotation.
 * @param old_state Immutable material-point input state selected by the element.
 * @param new_state Optional material-point output state selected by the element.
 * @return Derivatives of the global tangent with respect to the three parameters.
 */
std::array<Mat6, 3> SolidSection::tangent_rotation_derivatives(
    const Vec3&                position_reference,
    const Mat3&                additional_rotation,
    const std::array<Mat3, 3>& additional_rotation_derivatives,
    const Precision*           old_state,
    Precision*                 new_state
) const {
    // Validate the linear constitutive formulation used by this sensitivity.
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");

    auto elasticity = material_->elasticity();

    // Evaluate the orientation-independent local material tangent at zero strain.
    VolumeStrainLinearized zero_strain;
    VolumeStressCauchy     zero_stress;
    Mat6                   tangent_material;

    elasticity->evaluate(
        zero_strain,
        old_state,
        new_state,
        zero_stress,
        &tangent_material
    );

    // Compose the section basis and additional rotation:
    //
    //     Q = Q_section R
    //     dQ/dq_i = Q_section dR/dq_i.
    const Mat3 section_basis  = section_orientation_basis(position_reference);
    const Mat3 material_basis = section_basis * additional_rotation;

    // Build the current strain and stress transformation operators.
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    // Differentiate T_strain column-by-column by transforming one unit global
    // engineering strain tensor through Q^T E Q.
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

    // Differentiate T_stress column-by-column from Q S Q^T for one unit local
    // physical stress tensor.
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

    // Apply the product rule independently to all three supplied rotation
    // directions.
    for (Index i = 0; i < 3; ++i) {
        const Mat3 material_basis_derivative =
            section_basis * additional_rotation_derivatives[i];
        const Mat6 strain_derivative =
            strain_transform_derivative(material_basis_derivative);
        const Mat6 stress_derivative =
            stress_transform_derivative(material_basis_derivative);

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
