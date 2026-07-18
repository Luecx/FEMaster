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
#include <utility>

namespace fem {

IntegratedShellSection::IntegratedShellSection(material::Material::Ptr    material,
                                               model::ElementRegion::Ptr  region,
                                               Precision                  thickness,
                                               cos::CoordinateSystem::Ptr orientation)
    : ShellSection(std::move(material), std::move(region), thickness, std::move(orientation)) {}

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
