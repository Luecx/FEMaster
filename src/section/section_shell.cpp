/**
 * @file section_shell.cpp
 * @brief Implements common shell section behavior and basis transformations.
 *
 * @see src/section/section_shell.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "section_shell.h"

#include "../core/logging.h"

#include <Eigen/Geometry>

#include <sstream>
#include <utility>

namespace fem {

ShellSection::ShellSection(material::Material::Ptr    material,
                           model::ElementRegion::Ptr  region,
                           Precision                  thickness,
                           cos::CoordinateSystem::Ptr orientation)
    : thickness_  (thickness),
      orientation_(std::move(orientation)) {
    logging::error(thickness_ > Precision(0), "ShellSection: thickness must be positive");

    material_ = std::move(material);
    region_   = std::move(region);
}

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
} // namespace fem
