/**
 * @file section_shell_abd.cpp
 * @brief Implements the shell section based on prescribed ABD and shear matrices.
 *
 * @see src/section/section_shell_abd.h
 * @author Finn Eggers
 * @date 18.07.2026
 */

#include "section_shell_abd.h"

#include "../core/logging.h"

#include <Eigen/Cholesky>

#include <sstream>
#include <utility>

namespace fem {

ABDShellSection::ABDShellSection(material::Material::Ptr    material,
                                 model::ElementRegion::Ptr  region,
                                 Precision                  thickness,
                                 const Mat6&                abd,
                                 const Mat2&                shear,
                                 cos::CoordinateSystem::Ptr orientation)
    : ShellSection(std::move(material), std::move(region), thickness, std::move(orientation)),
      abd_         (abd),
      shear_       (shear) {
    logging::error(abd_.allFinite(),
        "ABDShellSection: ABD matrix must contain only finite values");
    logging::error(shear_.allFinite(),
        "ABDShellSection: shear matrix must contain only finite values");
    logging::error(abd_.isApprox(abd_.transpose()),
        "ABDShellSection: ABD matrix must be symmetric");
    logging::error(shear_.isApprox(shear_.transpose()),
        "ABDShellSection: shear matrix must be symmetric");

    const Eigen::LLT<Mat6> abd_factorization(abd_);
    const Eigen::LLT<Mat2> shear_factorization(shear_);

    logging::error(abd_factorization.info() == Eigen::Success,
        "ABDShellSection: ABD matrix must be positive definite");
    logging::error(shear_factorization.info() == Eigen::Success,
        "ABDShellSection: shear matrix must be positive definite");
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
