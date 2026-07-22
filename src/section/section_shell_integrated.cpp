/**
 * @file section_shell_integrated.cpp
 * @brief Implements through-thickness shell material integration.
 *
 * @see src/section/section_shell_integrated.h
 * @author Finn Eggers
 * @date 18.07.2026
 */

#include "section_shell_integrated.h"

#include "../core/logging.h"
#include "../material/strain/shell_material_strain_green_lagrange.h"
#include "../material/strain/shell_material_strain_linearized.h"
#include "../material/stress/shell_material_stress_cauchy.h"
#include "../material/stress/shell_material_stress_pk2.h"

#include <array>
#include <cmath>
#include <utility>

namespace fem {

IntegratedShellSection::IntegratedShellSection(material::Material::Ptr    material,
                                               model::ElementRegion::Ptr  region,
                                               Precision                  thickness,
                                               cos::CoordinateSystem::Ptr orientation,
                                               Index                      csys_axis)
    : ShellSection(std::move(material),
                   std::move(region),
                   thickness,
                   std::move(orientation),
                   csys_axis) {}

/**
 * Evaluates the shell material model at one through-thickness point.
 *
 * The generalized shell strain is first expressed in the material basis when a
 * section orientation is present. Otherwise the element's local shell basis is
 * used for material evaluation and the final Cauchy stress is returned in
 * global components.
 *
 * Linearized material evaluation returns Cauchy stress directly. Green-Lagrange
 * evaluation returns PK2 stress, which is pushed forward with the supplied
 * deformation gradient before output.
 *
 * @param position_reference Physical reference position of the evaluation point.
 * @param shell_basis_global Orthonormal shell basis defining the input strain.
 * @param strain Generalized shell strain in the supplied shell basis.
 * @param z Physical through-thickness coordinate measured from the midsurface.
 * @param use_green_lagrange Select finite-strain material evaluation.
 * @param deformation_gradient Deformation gradient used for finite-strain recovery.
 * @return Cauchy stress in global components or in the projected material basis.
 */
VolumeStressCauchy IntegratedShellSection::compute_stress(const Vec3&                   position_reference,
                                                          const Mat3&                   shell_basis_global,
                                                          const ShellGeneralizedStrain& strain,
                                                          Precision                     z,
                                                          bool                          use_green_lagrange,
                                                          Mat3                          deformation_gradient) const {
    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;

    logging::error(material_ && material_->has_elasticity(),
        "IntegratedShellSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_shell_integration_green_lagrange() || !use_green_lagrange,
        "IntegratedShellSection material does not support Green-Lagrange shell material-point evaluation");
    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell material-point evaluation");
    logging::error(material_->elasticity()->state_size() == 0,
        "IntegratedShellSection does not yet provide material-point state storage");

    // Select the basis in which the material model receives strain components.
    Mat3                   stress_basis_global = shell_basis_global;
    ShellGeneralizedStrain strain_material     = strain;

    if (orientation_) {
        const Mat2 output_to_shell  = output_to_shell_rotation(position_reference, shell_basis_global);
        const Mat8 strain_transform = ShellGeneralizedStrain::transformation(output_to_shell);

        stress_basis_global = output_basis_global(position_reference, shell_basis_global);
        strain_material     = ShellGeneralizedStrain(strain_transform * strain.values());
    }

    // Reconstruct the material-point shell strain at the requested thickness
    // coordinate in the selected material evaluation basis.
    ShellMaterialStrain material_strain;
    material_strain[MaterialStrainComponent::XX] =
        strain_material[GeneralizedStrainComponent::EpsilonXX]
        + z * strain_material[GeneralizedStrainComponent::KappaXX];
    material_strain[MaterialStrainComponent::YY] =
        strain_material[GeneralizedStrainComponent::EpsilonYY]
        + z * strain_material[GeneralizedStrainComponent::KappaYY];
    material_strain[MaterialStrainComponent::GammaXY] =
        strain_material[GeneralizedStrainComponent::GammaXY]
        + z * strain_material[GeneralizedStrainComponent::KappaXY];
    material_strain[MaterialStrainComponent::GammaXZ] = strain_material[GeneralizedStrainComponent::GammaXZ];
    material_strain[MaterialStrainComponent::GammaYZ] = strain_material[GeneralizedStrainComponent::GammaYZ];

    // Convert the five plane-stress shell components into a full symmetric
    // tensor in the selected material evaluation basis.
    auto shell_stress_tensor = [](const ShellMaterialStress& stress) {
        Mat3 tensor = Mat3::Zero();
        tensor(0, 0) = stress[ShellMaterialStress::Component::XX];
        tensor(1, 1) = stress[ShellMaterialStress::Component::YY];
        tensor(0, 1) = stress[ShellMaterialStress::Component::XY];
        tensor(1, 0) = stress[ShellMaterialStress::Component::XY];
        tensor(0, 2) = stress[ShellMaterialStress::Component::XZ];
        tensor(2, 0) = stress[ShellMaterialStress::Component::XZ];
        tensor(1, 2) = stress[ShellMaterialStress::Component::YZ];
        tensor(2, 1) = stress[ShellMaterialStress::Component::YZ];
        return tensor;
    };

    Mat5 material_tangent;

    if (!use_green_lagrange) {
        ShellMaterialStrainLinearized material_strain_linearized(material_strain.values());
        ShellMaterialStressCauchy     material_stress_cauchy;

        // Linearized shell material stresses are already Cauchy stresses in the
        // selected material evaluation basis.
        material_->elasticity()->evaluate(
            material_strain_linearized,
            nullptr,
            nullptr,
            material_stress_cauchy,
            material_tangent
        );

        const Mat3 cauchy_tensor = shell_stress_tensor(material_stress_cauchy);
        if (orientation_) {
            return VolumeStressCauchy(cauchy_tensor);
        }

        return VolumeStressCauchy(stress_basis_global * cauchy_tensor * stress_basis_global.transpose());
    }

    ShellMaterialStrainGreenLagrange material_strain_gl(material_strain.values());
    ShellMaterialStressPK2           material_stress_pk2;

    // Finite-strain material stresses are PK2 stresses and must be pushed
    // forward before Cauchy stress output.
    material_->elasticity()->evaluate(
        material_strain_gl,
        nullptr,
        nullptr,
        material_stress_pk2,
        material_tangent
    );

    const Precision J = deformation_gradient.determinant();
    logging::error(J > Precision(0) && std::isfinite(J),
        "IntegratedShellSection: invalid deformation gradient during stress recovery, J = ", J);

    const Mat3 second_pk_basis  = shell_stress_tensor(material_stress_pk2);
    const Mat3 second_pk_global = stress_basis_global * second_pk_basis * stress_basis_global.transpose();
    const Mat3 cauchy_global =
        (deformation_gradient * second_pk_global * deformation_gradient.transpose()) / J;

    if (!orientation_) {
        return VolumeStressCauchy(cauchy_global);
    }

    return VolumeStressCauchy(cauchy_global).transformed(Mat3::Identity(), stress_basis_global);
}

void IntegratedShellSection::evaluate_material(const ShellGeneralizedStrain& strain,
                                               bool                          use_green_lagrange,
                                               ShellStressResultants&        resultants,
                                               Mat8&                         tangent) const {

    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;
    using ResultantComponent         = ShellStressResultants::Component;

    logging::error(material_ && material_->has_elasticity(),
        "IntegratedShellSection requires a material with elasticity");
    logging::error(material_->elasticity()->supports_shell_integration_green_lagrange() || !use_green_lagrange,
        "IntegratedShellSection material does not support Green-Lagrange shell material-point evaluation");
    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell material-point evaluation");
    logging::error(material_->elasticity()->state_size() == 0,
        "IntegratedShellSection does not yet provide material-point state storage");

    const Precision shear_correction = Precision(5) / Precision(6);

    // 5 point simpson rule
    const std::array<Precision, 5> points {
        Precision(-1),
        Precision(-0.5),
        Precision(0),
        Precision(0.5),
        Precision(1)
    };
    const std::array<Precision, 5> weights {
        Precision(1) / Precision(6),
        Precision(4) / Precision(6),
        Precision(2) / Precision(6),
        Precision(4) / Precision(6),
        Precision(1) / Precision(6)
    };

    resultants.values().setZero();
    tangent.setZero();

    for (Index mp = 0; mp < 5; ++mp) {
        const Precision z = Precision(0.5) * thickness_ * points[mp];
        const Precision w = Precision(0.5) * thickness_ * weights[mp];

        ShellMaterialStrain material_strain;
        material_strain[MaterialStrainComponent::XX] =
            strain[GeneralizedStrainComponent::EpsilonXX]
            + z * strain[GeneralizedStrainComponent::KappaXX];
        material_strain[MaterialStrainComponent::YY] =
            strain[GeneralizedStrainComponent::EpsilonYY]
            + z * strain[GeneralizedStrainComponent::KappaYY];
        material_strain[MaterialStrainComponent::GammaXY] =
            strain[GeneralizedStrainComponent::GammaXY]
            + z * strain[GeneralizedStrainComponent::KappaXY];
        material_strain[MaterialStrainComponent::GammaXZ] = strain[GeneralizedStrainComponent::GammaXZ];
        material_strain[MaterialStrainComponent::GammaYZ] = strain[GeneralizedStrainComponent::GammaYZ];

        ShellMaterialStress material_stress;
        Mat5                material_tangent;

        if (use_green_lagrange) {
            ShellMaterialStrainGreenLagrange material_strain_gl(material_strain.values());
            ShellMaterialStressPK2           material_stress_pk2;

            material_->elasticity()->evaluate(
                material_strain_gl,
                nullptr,
                nullptr,
                material_stress_pk2,
                material_tangent
            );

            material_stress.values() = material_stress_pk2.values();
        } else {
            ShellMaterialStrainLinearized material_strain_linearized(material_strain.values());
            ShellMaterialStressCauchy     material_stress_cauchy;

            material_->elasticity()->evaluate(
                material_strain_linearized,
                nullptr,
                nullptr,
                material_stress_cauchy,
                material_tangent
            );

            material_stress.values() = material_stress_cauchy.values();
        }

        resultants[ResultantComponent::NXX] += w                    * material_stress[MaterialStressComponent::XX];
        resultants[ResultantComponent::NYY] += w                    * material_stress[MaterialStressComponent::YY];
        resultants[ResultantComponent::NXY] += w                    * material_stress[MaterialStressComponent::XY];
        resultants[ResultantComponent::MXX] += w * z                * material_stress[MaterialStressComponent::XX];
        resultants[ResultantComponent::MYY] += w * z                * material_stress[MaterialStressComponent::YY];
        resultants[ResultantComponent::MXY] += w * z                * material_stress[MaterialStressComponent::XY];
        resultants[ResultantComponent::QX]  += w * shear_correction * material_stress[MaterialStressComponent::XZ];
        resultants[ResultantComponent::QY]  += w * shear_correction * material_stress[MaterialStressComponent::YZ];

        StaticMatrix<5, 8> strain_map = StaticMatrix<5, 8>::Zero();
        strain_map(static_cast<Index>(MaterialStrainComponent::XX),
                   static_cast<Index>(GeneralizedStrainComponent::EpsilonXX)) = Precision(1);
        strain_map(static_cast<Index>(MaterialStrainComponent::XX),
                   static_cast<Index>(GeneralizedStrainComponent::KappaXX)) = z;
        strain_map(static_cast<Index>(MaterialStrainComponent::YY),
                   static_cast<Index>(GeneralizedStrainComponent::EpsilonYY)) = Precision(1);
        strain_map(static_cast<Index>(MaterialStrainComponent::YY),
                   static_cast<Index>(GeneralizedStrainComponent::KappaYY)) = z;
        strain_map(static_cast<Index>(MaterialStrainComponent::GammaXY),
                   static_cast<Index>(GeneralizedStrainComponent::GammaXY)) = Precision(1);
        strain_map(static_cast<Index>(MaterialStrainComponent::GammaXY),
                   static_cast<Index>(GeneralizedStrainComponent::KappaXY)) = z;
        strain_map(static_cast<Index>(MaterialStrainComponent::GammaXZ),
                   static_cast<Index>(GeneralizedStrainComponent::GammaXZ)) = Precision(1);
        strain_map(static_cast<Index>(MaterialStrainComponent::GammaYZ),
                   static_cast<Index>(GeneralizedStrainComponent::GammaYZ)) = Precision(1);

        StaticMatrix<8, 5> resultant_map = StaticMatrix<8, 5>::Zero();
        resultant_map(static_cast<Index>(ResultantComponent::NXX),
                      static_cast<Index>(MaterialStressComponent::XX)) = Precision(1);
        resultant_map(static_cast<Index>(ResultantComponent::NYY),
                      static_cast<Index>(MaterialStressComponent::YY)) = Precision(1);
        resultant_map(static_cast<Index>(ResultantComponent::NXY),
                      static_cast<Index>(MaterialStressComponent::XY)) = Precision(1);
        resultant_map(static_cast<Index>(ResultantComponent::MXX),
                      static_cast<Index>(MaterialStressComponent::XX)) = z;
        resultant_map(static_cast<Index>(ResultantComponent::MYY),
                      static_cast<Index>(MaterialStressComponent::YY)) = z;
        resultant_map(static_cast<Index>(ResultantComponent::MXY),
                      static_cast<Index>(MaterialStressComponent::XY)) = z;
        resultant_map(static_cast<Index>(ResultantComponent::QX),
                      static_cast<Index>(MaterialStressComponent::XZ)) = shear_correction;
        resultant_map(static_cast<Index>(ResultantComponent::QY),
                      static_cast<Index>(MaterialStressComponent::YZ)) = shear_correction;

        tangent += w * resultant_map * material_tangent * strain_map;
    }
}
} // namespace fem
