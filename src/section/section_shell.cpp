/**
 * @file section_shell.cpp
 * @brief Implements shell section matrices and reporting.
 *
 * @see src/section/section_shell.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "section_shell.h"

#include "../core/logging.h"
#include "../material/strain/shell_material_strain_green_lagrange.h"
#include "../material/strain/shell_material_strain_linearized.h"
#include "../material/stress/shell_material_stress_cauchy.h"
#include "../material/stress/shell_material_stress_pk2.h"

#include <Eigen/Geometry>

#include <array>
#include <sstream>

namespace fem {

void ShellSection::evaluate(const Vec3&                   position_reference,
                            const Mat3&                   shell_basis_global,
                            const ShellGeneralizedStrain& strain,
                            bool                          use_green_lagrange,
                            ShellStressResultants&        resultants,
                            Mat8&                         tangent) const {
    if (!orientation_) {
        evaluate_material(strain, use_green_lagrange, resultants, tangent);
        return;
    }

    // Build the material basis in the supplied shell plane
    const Vec3 point_local      = orientation_->to_local(position_reference);
    const Mat3 orientation_axes = orientation_->get_axes(point_local);
    const Vec3 normal           = shell_basis_global.col(2).normalized();

    Vec3 material_e1 = orientation_axes.col(0);
    material_e1 -= normal * material_e1.dot(normal);

    if (material_e1.squaredNorm() < Precision(1e-24)) {
        material_e1 = shell_basis_global.col(0);
    }
    material_e1.normalize();

    const Vec3 material_e2 = normal.cross(material_e1).normalized();

    Mat3 material_basis_global;
    material_basis_global.col(0) = material_e1;
    material_basis_global.col(1) = material_e2;
    material_basis_global.col(2) = normal;

    // Transform generalized shell strains into the material basis
    Mat2 material_to_shell;
    material_to_shell(0, 0) = shell_basis_global.col(0).dot(material_basis_global.col(0));
    material_to_shell(1, 0) = shell_basis_global.col(1).dot(material_basis_global.col(0));
    material_to_shell(0, 1) = shell_basis_global.col(0).dot(material_basis_global.col(1));
    material_to_shell(1, 1) = shell_basis_global.col(1).dot(material_basis_global.col(1));

    const Precision a1 = material_to_shell(0, 0);
    const Precision a2 = material_to_shell(1, 0);
    const Precision b1 = material_to_shell(0, 1);
    const Precision b2 = material_to_shell(1, 1);

    Mat3 plane_strain_transform;
    plane_strain_transform <<
        a1 * a1,     a2 * a2,     a1 * a2,
        b1 * b1,     b2 * b2,     b1 * b2,
        2 * a1 * b1, 2 * a2 * b2, a1 * b2 + a2 * b1;

    using StrainComponent = ShellGeneralizedStrain::Component;

    const Index membrane_start = static_cast<Index>(StrainComponent::EpsilonXX);
    const Index curvature_start = static_cast<Index>(StrainComponent::KappaXX);
    const Index shear_start     = static_cast<Index>(StrainComponent::GammaXZ);

    Mat8 strain_transform = Mat8::Zero();
    strain_transform.template block<3, 3>(membrane_start, membrane_start) = plane_strain_transform;
    strain_transform.template block<3, 3>(curvature_start, curvature_start) = plane_strain_transform;
    strain_transform.template block<2, 2>(shear_start, shear_start) = material_to_shell.transpose();

    const Vec8                   strain_material_values = strain_transform * strain.values();
    const ShellGeneralizedStrain strain_material(strain_material_values);

    ShellStressResultants resultants_material;
    Mat8                  tangent_material;
    evaluate_material(
        strain_material,
        use_green_lagrange,
        resultants_material,
        tangent_material
    );

    // Return resultants and tangent in the supplied shell basis
    resultants.values() = strain_transform.transpose() * resultants_material.values();
    tangent             = strain_transform.transpose() * tangent_material * strain_transform;
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

void ABDShellSection::evaluate_material(const ShellGeneralizedStrain& strain,
                                        bool                          use_green_lagrange,
                                        ShellStressResultants&        resultants,
                                        Mat8&                         tangent) const {
    (void) use_green_lagrange;

    const Index membrane_row = static_cast<Index>(ShellStressResultants::Component::NXX);
    const Index membrane_col = static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);
    const Index shear_row    = static_cast<Index>(ShellStressResultants::Component::QX);
    const Index shear_col    = static_cast<Index>(ShellGeneralizedStrain::Component::GammaXZ);

    tangent.setZero();
    tangent.template block<6, 6>(membrane_row, membrane_col) = abd_;
    tangent.template block<2, 2>(shear_row, shear_col) = shear_;
    resultants.values() = tangent * strain.values();
}

void ShellSection::info() {
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material   : ", (material_    ? material_   ->name : "-"));
    logging::info(true, "   Region     : ", (region_      ? region_     ->name : "-"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
    logging::info(true, "   Thickness  : ", thickness_);
}

std::string ShellSection::str() const {
    std::ostringstream os;

    os << "ShellSection: t=" << thickness_
       << ", material="      << (material_     ? material_   ->name : std::string("-"))
       << ", orientation="   << (orientation_  ? orientation_->name : std::string("-"))
       << ", region="        << (region_       ? region_     ->name : std::string("-"))
       << " ("               << (region_       ? static_cast<int>(region_->size()) : 0) << ")";

    return os.str();
}

void ABDShellSection::info() {
    logging::info(true, "ABDShellSection:");
    logging::info(true, "   Material   : ", (material_    ? material_   ->name : "-"));
    logging::info(true, "   Region     : ", (region_      ? region_     ->name : "-"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
    logging::info(true, "   Thickness  : ", thickness_);
}

std::string ABDShellSection::str() const {
    std::ostringstream os;

    os << "ABDShellSection: t=" << thickness_
       << ", material="         << (material_    ? material_   ->name : std::string("-"))
       << ", orientation="      << (orientation_ ? orientation_->name : std::string("-"))
       << ", region="           << (region_      ? static_cast<int>(region_->size()) : 0);

    return os.str();
}
} // namespace fem
