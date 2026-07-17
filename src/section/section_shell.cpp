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

#include <array>
#include <sstream>

namespace fem {

void IntegratedShellSection::evaluate(const ShellGeneralizedStrain& strain,
                                      bool                          use_green_lagrange,
                                      ShellStressResultants&        resultants,
                                      Mat8&                         tangent) const {
    logging::error(material_ && material_->has_elasticity(),
                   "IntegratedShellSection requires a material with elasticity");
    if (use_green_lagrange) {
        logging::error(material_->elasticity()->supports_shell_integration_green_lagrange(),
                       "IntegratedShellSection material does not support Green-Lagrange shell material-point evaluation");
    } else {
        logging::error(material_->elasticity()->supports_shell_integration_linearized(),
                       "IntegratedShellSection material does not support linearized shell material-point evaluation");
    }
    logging::error(material_->elasticity()->state_size() == 0,
                   "IntegratedShellSection does not yet provide material-point state storage");

    const Precision shear_correction = Precision(5) / Precision(6);

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

        Vec5 material_strain_values;
        material_strain_values << strain.values()(0) + z * strain.values()(3),
                                  strain.values()(1) + z * strain.values()(4),
                                  strain.values()(2) + z * strain.values()(5),
                                  strain.values()(6),
                                  strain.values()(7);

        Vec5 material_stress_values;
        Mat5 material_tangent;

        if (use_green_lagrange) {
            ShellMaterialStrainGreenLagrange material_strain(material_strain_values);
            ShellMaterialStressPK2           material_stress;

            material_->elasticity()->evaluate(
                material_strain,
                nullptr,
                nullptr,
                material_stress,
                material_tangent
            );

            material_stress_values = material_stress.values();
        } else {
            ShellMaterialStrainLinearized material_strain(material_strain_values);
            ShellMaterialStressCauchy     material_stress;

            material_->elasticity()->evaluate(
                material_strain,
                nullptr,
                nullptr,
                material_stress,
                material_tangent
            );

            material_stress_values = material_stress.values();
        }

        resultants.values()(0) += w                    * material_stress_values(0);
        resultants.values()(1) += w                    * material_stress_values(1);
        resultants.values()(2) += w                    * material_stress_values(2);
        resultants.values()(3) += w * z                * material_stress_values(0);
        resultants.values()(4) += w * z                * material_stress_values(1);
        resultants.values()(5) += w * z                * material_stress_values(2);
        resultants.values()(6) += w * shear_correction * material_stress_values(3);
        resultants.values()(7) += w * shear_correction * material_stress_values(4);

        StaticMatrix<5, 8> strain_map = StaticMatrix<5, 8>::Zero();
        strain_map(0, 0) = Precision(1);
        strain_map(0, 3) = z;
        strain_map(1, 1) = Precision(1);
        strain_map(1, 4) = z;
        strain_map(2, 2) = Precision(1);
        strain_map(2, 5) = z;
        strain_map(3, 6) = Precision(1);
        strain_map(4, 7) = Precision(1);

        StaticMatrix<8, 5> resultant_map = StaticMatrix<8, 5>::Zero();
        resultant_map(0, 0) = Precision(1);
        resultant_map(1, 1) = Precision(1);
        resultant_map(2, 2) = Precision(1);
        resultant_map(3, 0) = z;
        resultant_map(4, 1) = z;
        resultant_map(5, 2) = z;
        resultant_map(6, 3) = shear_correction;
        resultant_map(7, 4) = shear_correction;

        tangent += w * resultant_map * material_tangent * strain_map;
    }
}

void ABDShellSection::evaluate(const ShellGeneralizedStrain& strain,
                               bool                          use_green_lagrange,
                               ShellStressResultants&        resultants,
                               Mat8&                         tangent) const {
    (void) use_green_lagrange;

    tangent.setZero();
    tangent.template block<6, 6>(0, 0) = abd_;
    tangent.template block<2, 2>(6, 6) = shear_;
    resultants.values() = tangent * strain.values();
}

bool ShellSection::has_density() const {
    return material_ && material_->has_density();
}

Precision ShellSection::get_density() const {
    if (has_density()) {
        return material_->get_density();
    }

    return Precision(0);
}

void ShellSection::info() {
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material   : ", (material_ ? material_->name : "-"));
    logging::info(true, "   Region     : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Thickness  : ", thickness_);
    logging::info(true, "   Density    : ", (has_density() ? std::to_string(get_density()) : "NO"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
}

std::string ShellSection::str() const {
    std::ostringstream os;

    os << "ShellSection: t=" << thickness_
       << ", material="      << (material_ ? material_->name : std::string("-"))
       << ", density="       << (has_density() ? std::to_string(get_density()) : std::string("-"))
       << ", orientation="   << (orientation_ ? orientation_->name : std::string("-"))
       << ", region="        << (region_ ? region_->name : std::string("-"))
       << " ("               << (region_ ? static_cast<int>(region_->size()) : 0) << ")";

    return os.str();
}

void ABDShellSection::info() {
    logging::info(true, "ABDShellSection:");
    logging::info(true, "   Material   : ", (material_ ? material_->name : "-"));
    logging::info(true, "   Region     : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Thickness  : ", thickness_);
    logging::info(true, "   Density    : ", (has_density() ? std::to_string(get_density()) : "NO"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
}

std::string ABDShellSection::str() const {
    std::ostringstream os;

    os << "ABDShellSection: t="
       << thickness_
       << ", material="
       << (material_ ? material_->name : std::string("-"))
       << ", density="
       << (has_density() ? std::to_string(get_density()) : std::string("-"))
       << ", orientation="
       << (orientation_ ? orientation_->name : std::string("-"))
       << ", region="
       << (region_ ? static_cast<int>(region_->size()) : 0);

    return os.str();
}
} // namespace fem
