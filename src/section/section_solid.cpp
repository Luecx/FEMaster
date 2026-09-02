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
 * The section requires a tangent for this overload and therefore passes its
 * address explicitly to the pointer-based material interface.
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
    // Validate the requested constitutive formulation before transformation.
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");

    // Resolve the elastic law assigned through the generic material definition.
    auto elasticity = material_->elasticity();

    // Compose the spatial section orientation with the element-provided
    // additional material rotation.
    const Mat3 material_basis =
        section_orientation_basis(position_reference) * additional_rotation;

    // Retain both engineering-Voigt transformation operators because the same
    // pair is required for the consistent global tangent.
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    // Materialize the Eigen expression as Vec6 before selecting the Voigt
    // constructor. Both VolumeStrainLinearized(Mat3) and (Vec6) are explicit,
    // so passing an unevaluated Eigen product leaves overload resolution
    // ambiguous with current Eigen versions.
    const Vec6 strain_material_values = strain_transform * strain_global.voigt();
    const VolumeStrainLinearized strain_material(strain_material_values);

    // Evaluate Cauchy stress and explicitly request the material tangent needed
    // by this section overload.
    VolumeStressCauchy stress_material;
    Mat6               tangent_material;

    elasticity->evaluate(
        strain_material,
        old_state,
        new_state,
        stress_material,
        &tangent_material
    );

    // Transform stress and the complete tangent back into global coordinates:
    //
    //     C_global = T_stress C_material T_strain.
    //
    // Again materialize the transformed Voigt vector explicitly to avoid an
    // ambiguous Mat3/Vec6 constructor selection from the Eigen product type.
    const Vec6 stress_global_values = stress_transform * stress_material.voigt();
    stress_global  = VolumeStressCauchy(stress_global_values);
    tangent_global = stress_transform * tangent_material * strain_transform;
}

/**
 * Evaluates Total-Lagrangian solid stress and tangent in global coordinates.
 *
 * This reference overload delegates to the optional-tangent implementation so
 * full Newton assembly and residual-only assembly share exactly one section
 * transformation path.
 *
 * @param position_reference Physical reference position of the material point.
 * @param additional_rotation Element-provided material rotation.
 * @param strain_global Green-Lagrange strain in global reference coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Material-point output state row.
 * @param stress_global PK2 stress in global reference coordinates.
 * @param tangent_global Consistent global tangent `dS/dE`.
 */
void SolidSection::evaluate(const Vec3&                      position_reference,
                            const Mat3&                      additional_rotation,
                            const VolumeStrainGreenLagrange& strain_global,
                            const Precision*                 old_state,
                            Precision*                       new_state,
                            VolumeStressPK2&                 stress_global,
                            Mat6&                            tangent_global) const {
    evaluate(
        position_reference,
        additional_rotation,
        strain_global,
        old_state,
        new_state,
        stress_global,
        &tangent_global
    );
}

/**
 * Evaluates Total-Lagrangian solid response with an optional tangent.
 *
 * Green-Lagrange strain is transformed from the global reference basis into the
 * optional material basis. PK2 stress is always evaluated and transformed back.
 * If `tangent_global` is null, the constitutive law receives the same null
 * tangent request and no tangent transformation is performed. This is the path
 * used by nonlinear residual and line-search evaluations.
 *
 * @param position_reference Physical reference position of the material point.
 * @param additional_rotation Element-provided rotation applied after the section basis.
 * @param strain_global Green-Lagrange strain in global reference coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Material-point output state row.
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
    // Validate the exact finite-strain constitutive pair before transformation.
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_green_lagrange(),
        "SolidSection material does not support Green-Lagrange volume evaluation");

    auto elasticity = material_->elasticity();

    // Build reference material basis and both engineering-Voigt transforms.
    const Mat3 material_basis =
        section_orientation_basis(position_reference) * additional_rotation;
    const Mat6 strain_transform =
        VolumeStrain::get_transformation_matrix(Mat3::Identity(), material_basis);
    const Mat6 stress_transform =
        VolumeStress::get_transformation_matrix(material_basis, Mat3::Identity());

    // Materialize the transformed strain vector before constructing the typed
    // Green-Lagrange wrapper. This selects the Vec6 constructor unambiguously.
    const Vec6 strain_material_values = strain_transform * strain_global.voigt();
    const VolumeStrainGreenLagrange strain_material(strain_material_values);

    // Forward tangent optionality into the material evaluation. The local tangent
    // object is touched only when global tangent output is requested.
    VolumeStressPK2 stress_material;
    Mat6            tangent_material;

    elasticity->evaluate(
        strain_material,
        old_state,
        new_state,
        stress_material,
        tangent_global != nullptr ? &tangent_material : nullptr
    );

    // Transform the mandatory PK2 stress back to the global reference basis.
    // The explicit Vec6 temporary avoids ambiguous conversion of an unevaluated
    // Eigen product to either the tensor or Voigt stress constructor.
    const Vec6 stress_global_values = stress_transform * stress_material.voigt();
    stress_global = VolumeStressPK2(stress_global_values);

    if (tangent_global != nullptr) {
        *tangent_global = stress_transform * tangent_material * strain_transform;
    }
}

/**
 * Computes the derivative of the transformed linear material tangent with
 * respect to three additional material-orientation parameters.
 *
 * The global tangent is
 *
 *     C_global = T_stress(Q) C_material T_strain(Q).
 *
 * Assuming the linear material operator itself is independent of the rotation
 * parameter, the product rule gives
 *
 *     dC_global/dq_i
 *         = dT_stress/dq_i C_material T_strain
 *         + T_stress C_material dT_strain/dq_i.
 *
 * The operation therefore explicitly requests the material tangent at zero
 * strain from the pointer-based constitutive interface.
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
    // Validate the linear elastic constitutive response used by this sensitivity.
    logging::error(material_ && material_->has_elasticity(),
        "SolidSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_volume_linearized(),
        "SolidSection material does not support linearized volume evaluation");

    auto elasticity = material_->elasticity();

    // Evaluate the local material tangent at zero strain. This sensitivity cannot
    // use a stress-only constitutive query because C_material is the differentiated
    // operator's central factor.
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

    // The material basis is composed of the prescribed section orientation and
    // the element-provided additional rotation:
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

    // Apply the product rule to each additional-rotation parameter.
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
