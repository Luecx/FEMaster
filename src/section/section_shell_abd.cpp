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

#include <cmath>
#include <sstream>
#include <utility>

namespace fem {

ABDShellSection::ABDShellSection(material::Material::Ptr    material,
                                 model::ElementRegion::Ptr  region,
                                 Precision                  thickness,
                                 const Mat6&                abd,
                                 const Mat2&                shear,
                                 cos::CoordinateSystem::Ptr orientation,
                                 Index                      csys_axis)
    : ShellSection(std::move(material),
                   std::move(region),
                   thickness,
                   std::move(orientation),
                   csys_axis),
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

/**
 * Reconstructs an equivalent through-thickness Cauchy stress from ABD
 * resultants.
 *
 * ABD sections do not store a layerwise material model. For stress output they
 * are therefore interpreted as one homogeneous equivalent shell layer: membrane
 * forces and bending moments form a linear stress distribution through the
 * thickness, and transverse shear resultants are distributed uniformly.
 *
 * Resultants are evaluated in the section output basis. Without an orientation
 * the reconstructed Cauchy tensor is returned in global components. With an
 * orientation the tensor is returned in the projected section material basis.
 *
 * @param position_reference Physical reference position of the evaluation point.
 * @param shell_basis_global Orthonormal shell basis defining the input strain.
 * @param strain Generalized shell strain in the supplied shell basis.
 * @param z Physical through-thickness coordinate measured from the midsurface.
 * @param use_green_lagrange Select finite-strain PK2-to-Cauchy recovery.
 * @param deformation_gradient Deformation gradient used for finite-strain recovery.
 * @return Equivalent Cauchy stress in the documented output basis.
 */
VolumeStressCauchy ABDShellSection::compute_stress(const Vec3&                   position_reference,
                                                   const Mat3&                   shell_basis_global,
                                                   const ShellGeneralizedStrain& strain,
                                                   Precision                     z,
                                                   bool                          use_green_lagrange,
                                                   Mat3                          deformation_gradient) const {
    const Precision h = thickness_;

    // Evaluate resultants in the same output basis used by shell section forces.
    const ShellStressResultants resultants = compute_resultants(
        position_reference,
        shell_basis_global,
        strain,
        use_green_lagrange
    );

    const Vec3 plane_stress =
        resultants.membrane() / h
        + z * (Precision(12) / (h * h * h)) * resultants.moments();
    const Vec2 shear_stress = resultants.transverse_shear() / h;

    VolumeStressCauchy stress_output;
    stress_output[VolumeStress::Component::XX] = plane_stress(0);
    stress_output[VolumeStress::Component::YY] = plane_stress(1);
    stress_output[VolumeStress::Component::ZZ] = Precision(0);
    stress_output[VolumeStress::Component::YZ] = shear_stress(1);
    stress_output[VolumeStress::Component::XZ] = shear_stress(0);
    stress_output[VolumeStress::Component::XY] = plane_stress(2);

    const Mat3 output_basis = output_basis_global(position_reference, shell_basis_global);

    if (!use_green_lagrange) {
        if (orientation_) {
            return stress_output;
        }

        const Mat3 cauchy_global = output_basis * stress_output.tensor() * output_basis.transpose();
        return VolumeStressCauchy(cauchy_global);
    }

    // Treat the equivalent ABD stress as a PK2 stress in the reference output
    // basis and push it forward for Cauchy stress output.
    const Precision J = deformation_gradient.determinant();
    logging::error(J > Precision(0) && std::isfinite(J),
        "ABDShellSection: invalid deformation gradient during stress recovery, J = ", J);

    const Mat3 second_pk_global = output_basis * stress_output.tensor() * output_basis.transpose();
    const Mat3 cauchy_global =
        (deformation_gradient * second_pk_global * deformation_gradient.transpose()) / J;

    if (!orientation_) {
        return VolumeStressCauchy(cauchy_global);
    }

    return VolumeStressCauchy(cauchy_global).transformed(Mat3::Identity(), output_basis);
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
    logging::info(orientation_ != nullptr, "   CSYSAXIS   : ", csys_axis_ + 1);
    logging::info(true, "   Thickness  : ", thickness_);
}

std::string ABDShellSection::str() const {
    std::ostringstream os;

    os << "ABDShellSection: t=" << thickness_
       << ", material="         << (material_    ? material_   ->name : std::string("-"))
       << ", orientation="      << (orientation_ ? orientation_->name : std::string("-"))
       << ", csysaxis="         << (orientation_ ? std::to_string(csys_axis_ + 1) : std::string("-"))
       << ", region="           << (region_      ? static_cast<int>(region_->size()) : 0);

    return os.str();
}
} // namespace fem
