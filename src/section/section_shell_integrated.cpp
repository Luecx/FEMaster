/**
 * @file section_shell_integrated.cpp
 * @brief Implements explicit through-thickness shell material integration.
 *
 * Generalized membrane strains, curvatures and transverse-shear strains are
 * reconstructed at five Simpson points. The material model returns stress and a
 * five-component consistent tangent at every point. These quantities are
 * integrated directly into generalized resultants and an eight-component
 * section tangent.
 *
 * @see IntegratedShellSection
 *
 * @author Finn Eggers
 * @date 22.07.2026
 */

#include "section_shell_integrated.h"

#include "../core/logging.h"
#include "../material/strain/shell_material_strain_green_lagrange.h"
#include "../material/strain/shell_material_strain_linearized.h"
#include "../material/stress/shell_material_stress_cauchy.h"
#include "../material/stress/shell_material_stress_pk2.h"

#include <Eigen/LU>

#include <array>
#include <cmath>
#include <utility>

namespace fem {

namespace {

// Composite Simpson coordinates on the normalized thickness interval [-1,1].
constexpr std::array<Precision, 5> simpson_points {
    Precision(-1),
    Precision(-0.5),
    Precision(0),
    Precision(0.5),
    Precision(1)
};

// Composite Simpson weights on [-1,1]. Multiplication by h/2 maps the rule to
// the physical interval [-h/2,h/2].
constexpr std::array<Precision, 5> simpson_weights {
    Precision(1) / Precision(6),
    Precision(4) / Precision(6),
    Precision(2) / Precision(6),
    Precision(4) / Precision(6),
    Precision(1) / Precision(6)
};

} // namespace

IntegratedShellSection::IntegratedShellSection(
    material::Material::Ptr    material,
    model::ElementRegion::Ptr  region,
    Precision                  thickness,
    cos::CoordinateSystem::Ptr orientation,
    Index                      csys_axis
)
    : ShellSection(
          std::move(material),
          std::move(region),
          thickness,
          std::move(orientation),
          csys_axis
      ) {}

void IntegratedShellSection::evaluate(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    bool                          use_green_lagrange,
    ShellStressResultants&        resultants_shell,
    Mat8&                         tangent_shell
) const {
    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;
    using ResultantComponent         = ShellStressResultants::Component;

    // Generalized shell strain components.
    constexpr GeneralizedStrainComponent EXX = GeneralizedStrainComponent::EpsilonXX;
    constexpr GeneralizedStrainComponent EYY = GeneralizedStrainComponent::EpsilonYY;
    constexpr GeneralizedStrainComponent GXY = GeneralizedStrainComponent::GammaXY;
    constexpr GeneralizedStrainComponent KXX = GeneralizedStrainComponent::KappaXX;
    constexpr GeneralizedStrainComponent KYY = GeneralizedStrainComponent::KappaYY;
    constexpr GeneralizedStrainComponent KXY = GeneralizedStrainComponent::KappaXY;
    constexpr GeneralizedStrainComponent GXZ = GeneralizedStrainComponent::GammaXZ;
    constexpr GeneralizedStrainComponent GYZ = GeneralizedStrainComponent::GammaYZ;

    // Material strain and stress components. Both use the same five-component
    // shell ordering, but are kept separate from the generalized shell indices.
    constexpr MaterialStrainComponent MXX = MaterialStrainComponent::XX;
    constexpr MaterialStrainComponent MYY = MaterialStrainComponent::YY;
    constexpr MaterialStrainComponent MXY = MaterialStrainComponent::GammaXY;
    constexpr MaterialStrainComponent MXZ = MaterialStrainComponent::GammaXZ;
    constexpr MaterialStrainComponent MYZ = MaterialStrainComponent::GammaYZ;

    constexpr MaterialStressComponent SXX = MaterialStressComponent::XX;
    constexpr MaterialStressComponent SYY = MaterialStressComponent::YY;
    constexpr MaterialStressComponent SXY = MaterialStressComponent::XY;
    constexpr MaterialStressComponent SXZ = MaterialStressComponent::XZ;
    constexpr MaterialStressComponent SYZ = MaterialStressComponent::YZ;

    // Generalized stress-resultant components.
    constexpr ResultantComponent NXX = ResultantComponent::NXX;
    constexpr ResultantComponent NYY = ResultantComponent::NYY;
    constexpr ResultantComponent NXY = ResultantComponent::NXY;
    constexpr ResultantComponent BXX = ResultantComponent::MXX;
    constexpr ResultantComponent BYY = ResultantComponent::MYY;
    constexpr ResultantComponent BXY = ResultantComponent::MXY;
    constexpr ResultantComponent QX  = ResultantComponent::QX;
    constexpr ResultantComponent QY  = ResultantComponent::QY;

    // The integrated formulation requires an elasticity law supporting the
    // strain measure selected by the shell kinematics.
    logging::error(material_ && material_->has_elasticity(),
        "IntegratedShellSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_shell_integration_green_lagrange() || !use_green_lagrange,
        "IntegratedShellSection material does not support Green-Lagrange shell evaluation");
    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell evaluation");
    logging::error(material_->elasticity()->state_size() == 0,
        "IntegratedShellSection does not yet provide material-point state storage");

    // Use the geometric shell basis when no orientation exists. Otherwise use
    // the projected coordinate-system basis for all material evaluations.
    const Mat3 section_basis_global = orientation_
        ? stress_basis(position_reference, shell_basis_global)
        : shell_basis_global;

    // Express section axes in shell coordinates and rotate generalized strains
    // into the material basis.
    const Mat2 section_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * section_basis_global.template block<3, 2>(0, 0);

    const Mat8 strain_shell_to_section =
        ShellGeneralizedStrain::transformation(section_axes_in_shell);

    const ShellGeneralizedStrain strain_section(
        strain_shell_to_section * strain_shell.values()
    );

    // Standard Reissner-Mindlin correction for the assumed constant transverse
    // shear strain through the thickness.
    const Precision shear_correction = Precision(5) / Precision(6);

    ShellStressResultants resultants_section;
    Mat8                  tangent_section;

    resultants_section.values().setZero();
    tangent_section.setZero();

    // Integrate all material points through the physical thickness.
    for (Index mp = 0; mp < 5; ++mp) {
        // Map normalized coordinate and weight to [-h/2,h/2].
        const Precision z = Precision(0.5) * thickness_ * simpson_points[mp];
        const Precision w = Precision(0.5) * thickness_ * simpson_weights[mp];

        // Reconstruct material-point strains:
        //
        //     epsilon(z) = epsilon_0 + z kappa.
        //
        // Transverse shear remains constant in the present shell kinematics.
        ShellMaterialStrain material_strain;
        material_strain[MXX] = strain_section[EXX] + z * strain_section[KXX];
        material_strain[MYY] = strain_section[EYY] + z * strain_section[KYY];
        material_strain[MXY] = strain_section[GXY] + z * strain_section[KXY];
        material_strain[MXZ] = strain_section[GXZ];
        material_strain[MYZ] = strain_section[GYZ];

        ShellMaterialStress material_stress;
        Mat5                material_tangent;

        // Evaluate either PK2 stress conjugate to Green-Lagrange strain or
        // Cauchy stress for the linearized shell response.
        if (use_green_lagrange) {
            const ShellMaterialStrainGreenLagrange material_strain_gl(material_strain.values());
            ShellMaterialStressPK2                 material_stress_pk2;

            material_->elasticity()->evaluate(
                material_strain_gl,
                nullptr,
                nullptr,
                material_stress_pk2,
                material_tangent
            );

            material_stress.values() = material_stress_pk2.values();
        } else {
            const ShellMaterialStrainLinearized material_strain_linearized(material_strain.values());
            ShellMaterialStressCauchy           material_stress_cauchy;

            material_->elasticity()->evaluate(
                material_strain_linearized,
                nullptr,
                nullptr,
                material_stress_cauchy,
                material_tangent
            );

            material_stress.values() = material_stress_cauchy.values();
        }

        // Integrate membrane forces.
        resultants_section[NXX] += w * material_stress[SXX];
        resultants_section[NYY] += w * material_stress[SYY];
        resultants_section[NXY] += w * material_stress[SXY];

        // Integrate bending moments with the physical thickness coordinate z.
        resultants_section[BXX] += w * z * material_stress[SXX];
        resultants_section[BYY] += w * z * material_stress[SYY];
        resultants_section[BXY] += w * z * material_stress[SXY];

        // Integrate corrected transverse shear forces.
        resultants_section[QX] += w * shear_correction * material_stress[SXZ];
        resultants_section[QY] += w * shear_correction * material_stress[SYZ];

        // Derivative of material-point strain with respect to generalized shell
        // strain at the current thickness coordinate.
        StaticMatrix<5, 8> strain_map = StaticMatrix<5, 8>::Zero();

        strain_map((Index)MXX, (Index)EXX) = Precision(1);
        strain_map((Index)MXX, (Index)KXX) = z;
        strain_map((Index)MYY, (Index)EYY) = Precision(1);
        strain_map((Index)MYY, (Index)KYY) = z;
        strain_map((Index)MXY, (Index)GXY) = Precision(1);
        strain_map((Index)MXY, (Index)KXY) = z;
        strain_map((Index)MXZ, (Index)GXZ) = Precision(1);
        strain_map((Index)MYZ, (Index)GYZ) = Precision(1);

        // Derivative of generalized resultants with respect to material stress.
        StaticMatrix<8, 5> resultant_map = StaticMatrix<8, 5>::Zero();

        resultant_map((Index)NXX, (Index)SXX) = Precision(1);
        resultant_map((Index)NYY, (Index)SYY) = Precision(1);
        resultant_map((Index)NXY, (Index)SXY) = Precision(1);
        resultant_map((Index)BXX, (Index)SXX) = z;
        resultant_map((Index)BYY, (Index)SYY) = z;
        resultant_map((Index)BXY, (Index)SXY) = z;
        resultant_map((Index)QX,  (Index)SXZ) = shear_correction;
        resultant_map((Index)QY,  (Index)SYZ) = shear_correction;

        // Chain-rule integration of the generalized consistent tangent:
        //
        //     H += w R_sigma C_material R_epsilon.
        tangent_section.noalias() += w * resultant_map * material_tangent * strain_map;
    }

    // With no orientation, section and shell bases are identical.
    if (!orientation_) {
        resultants_shell = resultants_section;
        tangent_shell    = tangent_section;
        return;
    }

    // Rotate the integrated physical resultants back to the geometric shell
    // basis. Resultants use a stress-type rather than engineering-strain
    // transformation.
    const Mat2 shell_axes_in_section = section_axes_in_shell.transpose();
    const Mat8 resultants_section_to_shell =
        ShellStressResultants::transformation(shell_axes_in_section);

    resultants_shell = ShellStressResultants(
        resultants_section_to_shell * resultants_section.values()
    );

    tangent_shell =
        resultants_section_to_shell
        * tangent_section
        * strain_shell_to_section;
}

VolumeStressCauchy IntegratedShellSection::evaluate_output_stress(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    Precision                     z,
    bool                          use_green_lagrange,
    const Mat3&                   deformation_gradient
) const {
    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;

    // Generalized shell strain components.
    constexpr GeneralizedStrainComponent EXX = GeneralizedStrainComponent::EpsilonXX;
    constexpr GeneralizedStrainComponent EYY = GeneralizedStrainComponent::EpsilonYY;
    constexpr GeneralizedStrainComponent GXY = GeneralizedStrainComponent::GammaXY;
    constexpr GeneralizedStrainComponent KXX = GeneralizedStrainComponent::KappaXX;
    constexpr GeneralizedStrainComponent KYY = GeneralizedStrainComponent::KappaYY;
    constexpr GeneralizedStrainComponent KXY = GeneralizedStrainComponent::KappaXY;
    constexpr GeneralizedStrainComponent GXZ = GeneralizedStrainComponent::GammaXZ;
    constexpr GeneralizedStrainComponent GYZ = GeneralizedStrainComponent::GammaYZ;

    // Material strain and stress components.
    constexpr MaterialStrainComponent MXX = MaterialStrainComponent::XX;
    constexpr MaterialStrainComponent MYY = MaterialStrainComponent::YY;
    constexpr MaterialStrainComponent MXY = MaterialStrainComponent::GammaXY;
    constexpr MaterialStrainComponent MXZ = MaterialStrainComponent::GammaXZ;
    constexpr MaterialStrainComponent MYZ = MaterialStrainComponent::GammaYZ;

    constexpr MaterialStressComponent SXX = MaterialStressComponent::XX;
    constexpr MaterialStressComponent SYY = MaterialStressComponent::YY;
    constexpr MaterialStressComponent SXY = MaterialStressComponent::XY;
    constexpr MaterialStressComponent SXZ = MaterialStressComponent::XZ;
    constexpr MaterialStressComponent SYZ = MaterialStressComponent::YZ;

    // Apply the same admissibility checks as the integrated generalized response.
    logging::error(material_ && material_->has_elasticity(),
        "IntegratedShellSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_shell_integration_green_lagrange() || !use_green_lagrange,
        "IntegratedShellSection material does not support Green-Lagrange shell evaluation");
    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell evaluation");
    logging::error(material_->elasticity()->state_size() == 0,
        "IntegratedShellSection does not yet provide material-point state storage");

    // Stress output is global without an orientation and section-local with one.
    const Mat3 output_basis_global = stress_basis(position_reference, shell_basis_global);

    // Material recovery always requires a tangential basis. Without an
    // orientation, use the geometric shell basis and rotate only the final tensor
    // to global components.
    const Mat3 recovery_basis_global = orientation_
        ? output_basis_global
        : shell_basis_global;

    // Rotate generalized strain into the material recovery basis.
    const Mat2 recovery_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * recovery_basis_global.template block<3, 2>(0, 0);

    const ShellGeneralizedStrain strain_material =
        strain_shell.transformed(recovery_axes_in_shell);

    // Reconstruct the five material-point shell strain components.
    ShellMaterialStrain material_strain;
    material_strain[MXX] = strain_material[EXX] + z * strain_material[KXX];
    material_strain[MYY] = strain_material[EYY] + z * strain_material[KYY];
    material_strain[MXY] = strain_material[GXY] + z * strain_material[KXY];
    material_strain[MXZ] = strain_material[GXZ];
    material_strain[MYZ] = strain_material[GYZ];

    // Convert five plane-stress shell components into a symmetric 3D tensor in
    // the recovery basis.
    auto shell_stress_tensor = [=](const ShellMaterialStress& stress) {
        Mat3 tensor = Mat3::Zero();

        tensor(0, 0) = stress[SXX];
        tensor(1, 1) = stress[SYY];
        tensor(0, 1) = stress[SXY];
        tensor(1, 0) = stress[SXY];
        tensor(0, 2) = stress[SXZ];
        tensor(2, 0) = stress[SXZ];
        tensor(1, 2) = stress[SYZ];
        tensor(2, 1) = stress[SYZ];

        return tensor;
    };

    // The material API currently requires a tangent output even though this
    // pointwise recovery path consumes only the stress.
    Mat5 material_tangent;

    if (!use_green_lagrange) {
        const ShellMaterialStrainLinearized material_strain_linearized(material_strain.values());
        ShellMaterialStressCauchy           material_stress_cauchy;

        material_->elasticity()->evaluate(
            material_strain_linearized,
            nullptr,
            nullptr,
            material_stress_cauchy,
            material_tangent
        );

        // Linearized material stress is already Cauchy stress in the recovery
        // basis. Transform it directly into the configured output basis.
        return VolumeStressCauchy(shell_stress_tensor(material_stress_cauchy))
            .transformed(recovery_basis_global, output_basis_global);
    }

    const ShellMaterialStrainGreenLagrange material_strain_gl(material_strain.values());
    ShellMaterialStressPK2                 material_stress_pk2;

    material_->elasticity()->evaluate(
        material_strain_gl,
        nullptr,
        nullptr,
        material_stress_pk2,
        material_tangent
    );

    // Validate the finite-strain push-forward before division by J.
    const Precision J = deformation_gradient.determinant();
    logging::error(
        J > Precision(0) && std::isfinite(J),
        "IntegratedShellSection: invalid deformation gradient during stress recovery, J = ",
        J
    );

    // Rotate PK2 stress from the reference material basis to global coordinates
    // and apply sigma = J^-1 F S F^T.
    const Mat3 second_pk_recovery = shell_stress_tensor(material_stress_pk2);
    const Mat3 second_pk_global   = recovery_basis_global * second_pk_recovery * recovery_basis_global.transpose();
    const Mat3 cauchy_global      = deformation_gradient * second_pk_global * deformation_gradient.transpose() / J;

    // Return global components without an orientation and section components
    // with an orientation.
    return VolumeStressCauchy(cauchy_global).transformed(Mat3::Identity(), output_basis_global);
}

} // namespace fem
