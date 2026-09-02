/**
 * @file section_shell_integrated.cpp
 * @brief Implements explicit through-thickness shell material integration.
 *
 * Generalized membrane strains, curvatures and transverse-shear strains are
 * reconstructed at five Simpson points. The material model returns stress and,
 * for generalized section assembly, a five-component consistent tangent at every
 * point. These quantities are integrated directly into generalized resultants
 * and an eight-component section tangent.
 *
 * A null output-state pointer denotes a state-neutral constitutive query. The
 * integrated section propagates that null target through every through-thickness
 * material point without performing pointer arithmetic on it. Stress-recovery
 * calls also pass a null tangent pointer because no material derivative is needed
 * for output conversion.
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

// Composite Simpson coordinates on the normalized thickness interval [-1,1]
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

/**
 * Constructs a five-point through-thickness integrated shell section.
 *
 * The base class validates thickness and orientation data. Constitutive support
 * for the requested strain measure is checked at evaluation time because the
 * owning material may provide only selected shell formulations.
 *
 * @param material Material evaluated at every through-thickness point.
 * @param region Element region receiving the section.
 * @param thickness Positive physical shell thickness.
 * @param orientation Optional coordinate system defining the material basis.
 * @param csys_axis Zero-based coordinate-system axis projected into the shell plane.
 */
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

/**
 * Integrates one shell integration point through the five physical material
 * points of the section.
 *
 * The state pointers address the first old/new material-point rows belonging to
 * the shell integration point. Consecutive rows are separated by
 * `material_state_stride` scalar components. `new_material_state == nullptr`
 * marks a state-neutral query and is forwarded as `nullptr` to every material
 * point instead of creating or advancing through-thickness trial history.
 *
 * Generalized membrane strains and curvatures reconstruct the in-plane material
 * strain through
 *
 *     epsilon(z) = epsilon_0 + z kappa,
 *
 * while transverse shear is constant through the thickness. Each material point
 * returns stress and an explicitly requested five-component material tangent.
 * Simpson integration then assembles membrane forces, bending moments,
 * transverse shear forces and the consistent generalized tangent. Oriented
 * quantities are transformed back into the geometric shell basis used by the
 * element.
 *
 * @param position_reference Physical reference position of the shell point.
 * @param shell_basis_global Geometric shell basis in global coordinates.
 * @param strain_shell Generalized strain in the geometric shell basis.
 * @param old_material_state First of the five through-thickness input state rows.
 * @param new_material_state Optional first through-thickness output state row.
 * @param material_state_stride Scalar distance between consecutive state rows.
 * @param use_green_lagrange Select PK2 or linearized Cauchy material evaluation.
 * @param resultants_shell Integrated resultants in the geometric shell basis.
 * @param tangent_shell Consistent generalized tangent in the geometric shell basis.
 */
void IntegratedShellSection::evaluate(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    const Precision*              old_material_state,
    Precision*                    new_material_state,
    Index                         material_state_stride,
    bool                          use_green_lagrange,
    ShellStressResultants&        resultants_shell,
    Mat8&                         tangent_shell
) const {
    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;
    using ResultantComponent         = ShellStressResultants::Component;

    constexpr GeneralizedStrainComponent EXX = GeneralizedStrainComponent::EpsilonXX;
    constexpr GeneralizedStrainComponent EYY = GeneralizedStrainComponent::EpsilonYY;
    constexpr GeneralizedStrainComponent GXY = GeneralizedStrainComponent::GammaXY;
    constexpr GeneralizedStrainComponent KXX = GeneralizedStrainComponent::KappaXX;
    constexpr GeneralizedStrainComponent KYY = GeneralizedStrainComponent::KappaYY;
    constexpr GeneralizedStrainComponent KXY = GeneralizedStrainComponent::KappaXY;
    constexpr GeneralizedStrainComponent GXZ = GeneralizedStrainComponent::GammaXZ;
    constexpr GeneralizedStrainComponent GYZ = GeneralizedStrainComponent::GammaYZ;

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

    constexpr ResultantComponent NXX = ResultantComponent::NXX;
    constexpr ResultantComponent NYY = ResultantComponent::NYY;
    constexpr ResultantComponent NXY = ResultantComponent::NXY;
    constexpr ResultantComponent BXX = ResultantComponent::MXX;
    constexpr ResultantComponent BYY = ResultantComponent::MYY;
    constexpr ResultantComponent BXY = ResultantComponent::MXY;
    constexpr ResultantComponent QX  = ResultantComponent::QX;
    constexpr ResultantComponent QY  = ResultantComponent::QY;

    // Validate the material and exact shell strain measure before constructing
    // any through-thickness constitutive state.
    logging::error(material_ && material_->has_elasticity(),
        "IntegratedShellSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_shell_integration_green_lagrange() || !use_green_lagrange,
        "IntegratedShellSection material does not support Green-Lagrange shell evaluation");
    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell evaluation");

    // Determine the material section basis. Without a prescribed orientation the
    // geometric shell basis itself is the material basis.
    const Mat3 section_basis_global = orientation_
        ? stress_basis(position_reference, shell_basis_global)
        : shell_basis_global;

    const Mat2 section_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * section_basis_global.template block<3, 2>(0, 0);

    // Transform generalized shell strain into material-section coordinates once
    // before integrating through the thickness.
    const Mat8 strain_shell_to_section =
        ShellGeneralizedStrain::transformation(section_axes_in_shell);
    const ShellGeneralizedStrain strain_section(
        strain_shell_to_section * strain_shell.values()
    );

    const Precision shear_correction = Precision(5) / Precision(6);

    ShellStressResultants resultants_section;
    Mat8                  tangent_section;

    resultants_section.values().setZero();
    tangent_section.setZero();

    // Integrate all five physical material points through the shell thickness.
    for (Index mp = 0; mp < 5; ++mp) {
        const Precision z = Precision(0.5) * thickness_ * simpson_points[mp];
        const Precision w = Precision(0.5) * thickness_ * simpson_weights[mp];

        // Reconstruct the five local material strain components at the current
        // thickness coordinate.
        ShellMaterialStrain material_strain;
        material_strain[MXX] = strain_section[EXX] + z * strain_section[KXX];
        material_strain[MYY] = strain_section[EYY] + z * strain_section[KYY];
        material_strain[MXY] = strain_section[GXY] + z * strain_section[KXY];
        material_strain[MXZ] = strain_section[GXZ];
        material_strain[MYZ] = strain_section[GYZ];

        // Resolve the material-point history rows without advancing a null output
        // pointer during state-neutral evaluations.
        const Precision* old_state = old_material_state + mp * material_state_stride;
        Precision* new_state = new_material_state
            ? new_material_state + mp * material_state_stride
            : nullptr;

        ShellMaterialStress material_stress;
        Mat5                material_tangent;

        // Generalized section stiffness requires the material derivative at every
        // Simpson point, so both branches pass its address explicitly.
        if (use_green_lagrange) {
            const ShellMaterialStrainGreenLagrange material_strain_gl(material_strain.values());
            ShellMaterialStressPK2                 material_stress_pk2;

            material_->elasticity()->evaluate(
                material_strain_gl,
                old_state,
                new_state,
                material_stress_pk2,
                &material_tangent
            );

            material_stress.values() = material_stress_pk2.values();
        } else {
            const ShellMaterialStrainLinearized material_strain_linearized(material_strain.values());
            ShellMaterialStressCauchy           material_stress_cauchy;

            material_->elasticity()->evaluate(
                material_strain_linearized,
                old_state,
                new_state,
                material_stress_cauchy,
                &material_tangent
            );

            material_stress.values() = material_stress_cauchy.values();
        }

        // Integrate membrane forces and bending moments from the in-plane stress
        // components. Bending resultants carry one additional factor z.
        resultants_section[NXX] += w * material_stress[SXX];
        resultants_section[NYY] += w * material_stress[SYY];
        resultants_section[NXY] += w * material_stress[SXY];

        resultants_section[BXX] += w * z * material_stress[SXX];
        resultants_section[BYY] += w * z * material_stress[SYY];
        resultants_section[BXY] += w * z * material_stress[SXY];

        // Apply the Reissner-Mindlin shear correction to transverse shear forces.
        resultants_section[QX] += w * shear_correction * material_stress[SXZ];
        resultants_section[QY] += w * shear_correction * material_stress[SYZ];

        // Map generalized strain increments to local material strain increments:
        //
        //     d epsilon_material = B_z d epsilon_generalized.
        StaticMatrix<5, 8> strain_map = StaticMatrix<5, 8>::Zero();
        strain_map((Index)MXX, (Index)EXX) = Precision(1);
        strain_map((Index)MXX, (Index)KXX) = z;
        strain_map((Index)MYY, (Index)EYY) = Precision(1);
        strain_map((Index)MYY, (Index)KYY) = z;
        strain_map((Index)MXY, (Index)GXY) = Precision(1);
        strain_map((Index)MXY, (Index)KXY) = z;
        strain_map((Index)MXZ, (Index)GXZ) = Precision(1);
        strain_map((Index)MYZ, (Index)GYZ) = Precision(1);

        // Map material stress increments into generalized resultant increments.
        StaticMatrix<8, 5> resultant_map = StaticMatrix<8, 5>::Zero();
        resultant_map((Index)NXX, (Index)SXX) = Precision(1);
        resultant_map((Index)NYY, (Index)SYY) = Precision(1);
        resultant_map((Index)NXY, (Index)SXY) = Precision(1);
        resultant_map((Index)BXX, (Index)SXX) = z;
        resultant_map((Index)BYY, (Index)SYY) = z;
        resultant_map((Index)BXY, (Index)SXY) = z;
        resultant_map((Index)QX,  (Index)SXZ) = shear_correction;
        resultant_map((Index)QY,  (Index)SYZ) = shear_correction;

        // Integrate the complete consistent generalized tangent contribution:
        //
        //     dR/dE_gen += w A_z C_material B_z.
        tangent_section.noalias() += w * resultant_map * material_tangent * strain_map;
    }

    // Without material orientation, section and shell bases are identical and no
    // result transformation is necessary.
    if (!orientation_) {
        resultants_shell = resultants_section;
        tangent_shell    = tangent_section;
        return;
    }

    // Transform integrated resultants and tangent back into geometric shell axes.
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

/**
 * Evaluates physical Cauchy stress at one requested thickness coordinate.
 *
 * The output point reuses the closest of the five Simpson material-point state
 * rows belonging to the in-plane shell integration point supplied by the
 * element. A null output-state pointer is propagated to the selected material
 * point so result recovery remains state-neutral.
 *
 * Material strain is reconstructed in the same section basis used during
 * generalized integration. Linearized constitutive output is already Cauchy
 * stress. Finite-strain output is evaluated as PK2 stress and pushed forward by
 * the supplied deformation gradient before conversion to the configured output
 * basis. No tangent participates in either recovery path, so both constitutive
 * calls explicitly pass `nullptr` for tangent output.
 *
 * @param position_reference Physical reference position of the shell point.
 * @param shell_basis_global Geometric shell basis in global coordinates.
 * @param strain_shell Generalized strain in the geometric shell basis.
 * @param old_material_state First of the five through-thickness input state rows.
 * @param new_material_state Optional first through-thickness output state row.
 * @param material_state_stride Scalar distance between consecutive state rows.
 * @param z Physical thickness coordinate measured from the midsurface.
 * @param use_green_lagrange Select PK2-to-Cauchy finite-strain recovery.
 * @param deformation_gradient Three-dimensional deformation gradient at `z`.
 * @return Physical Cauchy stress in the configured output basis.
 */
VolumeStressCauchy IntegratedShellSection::evaluate_output_stress(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    const Precision*              old_material_state,
    Precision*                    new_material_state,
    Index                         material_state_stride,
    Precision                     z,
    bool                          use_green_lagrange,
    const Mat3&                   deformation_gradient
) const {
    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;

    constexpr GeneralizedStrainComponent EXX = GeneralizedStrainComponent::EpsilonXX;
    constexpr GeneralizedStrainComponent EYY = GeneralizedStrainComponent::EpsilonYY;
    constexpr GeneralizedStrainComponent GXY = GeneralizedStrainComponent::GammaXY;
    constexpr GeneralizedStrainComponent KXX = GeneralizedStrainComponent::KappaXX;
    constexpr GeneralizedStrainComponent KYY = GeneralizedStrainComponent::KappaYY;
    constexpr GeneralizedStrainComponent KXY = GeneralizedStrainComponent::KappaXY;
    constexpr GeneralizedStrainComponent GXZ = GeneralizedStrainComponent::GammaXZ;
    constexpr GeneralizedStrainComponent GYZ = GeneralizedStrainComponent::GammaYZ;

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

    // Validate the material formulation selected for stress recovery.
    logging::error(material_ && material_->has_elasticity(),
        "IntegratedShellSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_shell_integration_green_lagrange() || !use_green_lagrange,
        "IntegratedShellSection material does not support Green-Lagrange shell evaluation");
    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell evaluation");

    // Determine both recovery and configured output bases. Without an explicit
    // orientation, material recovery is performed directly in shell coordinates.
    const Mat3 output_basis_global = stress_basis(position_reference, shell_basis_global);
    const Mat3 recovery_basis_global = orientation_
        ? output_basis_global
        : shell_basis_global;

    const Mat2 recovery_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * recovery_basis_global.template block<3, 2>(0, 0);

    const ShellGeneralizedStrain strain_material =
        strain_shell.transformed(recovery_axes_in_shell);

    // Reconstruct material strain at the requested physical thickness coordinate.
    ShellMaterialStrain material_strain;
    material_strain[MXX] = strain_material[EXX] + z * strain_material[KXX];
    material_strain[MYY] = strain_material[EYY] + z * strain_material[KYY];
    material_strain[MXY] = strain_material[GXY] + z * strain_material[KXY];
    material_strain[MXZ] = strain_material[GXZ];
    material_strain[MYZ] = strain_material[GYZ];

    // Select the nearest persistent Simpson material point because output
    // coordinates need not coincide with constitutive thickness points.
    Index     state_mp       = 0;
    Precision state_distance =
        std::abs(z - Precision(0.5) * thickness_ * simpson_points[0]);

    for (Index mp = 1; mp < 5; ++mp) {
        const Precision z_mp     = Precision(0.5) * thickness_ * simpson_points[mp];
        const Precision distance = std::abs(z - z_mp);

        if (distance < state_distance) {
            state_mp       = mp;
            state_distance = distance;
        }
    }

    const Precision* old_state = old_material_state + state_mp * material_state_stride;
    Precision* new_state = new_material_state
        ? new_material_state + state_mp * material_state_stride
        : nullptr;

    // Embed the five shell material stress components into a symmetric 3D tensor
    // for subsequent basis transformation and finite-strain push-forward.
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

    if (!use_green_lagrange) {
        // Linearized recovery already returns physical Cauchy stress. No material
        // tangent is required for result output.
        const ShellMaterialStrainLinearized material_strain_linearized(material_strain.values());
        ShellMaterialStressCauchy           material_stress_cauchy;

        material_->elasticity()->evaluate(
            material_strain_linearized,
            old_state,
            new_state,
            material_stress_cauchy,
            nullptr
        );

        return VolumeStressCauchy(shell_stress_tensor(material_stress_cauchy))
            .transformed(recovery_basis_global, output_basis_global);
    }

    // Finite-strain recovery obtains PK2 stress only and then pushes it forward
    // with the supplied deformation gradient.
    const ShellMaterialStrainGreenLagrange material_strain_gl(material_strain.values());
    ShellMaterialStressPK2                 material_stress_pk2;

    material_->elasticity()->evaluate(
        material_strain_gl,
        old_state,
        new_state,
        material_stress_pk2,
        nullptr
    );

    const Precision J = deformation_gradient.determinant();
    logging::error(J > Precision(0) && std::isfinite(J),
        "IntegratedShellSection: invalid deformation gradient during stress recovery, J = ", J);

    // Transform the recovered PK2 stress to global reference coordinates and
    // apply the standard push-forward
    //
    //     sigma = F S F^T / J.
    const Mat3 second_pk_recovery = shell_stress_tensor(material_stress_pk2);
    const Mat3 second_pk_global =
        recovery_basis_global * second_pk_recovery * recovery_basis_global.transpose();
    const Mat3 cauchy_global =
        deformation_gradient * second_pk_global * deformation_gradient.transpose() / J;

    return VolumeStressCauchy(cauchy_global)
        .transformed(Mat3::Identity(), output_basis_global);
}

} // namespace fem