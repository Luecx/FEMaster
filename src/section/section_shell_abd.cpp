/**
 * @file section_shell_abd.cpp
 * @brief Implements the shell section based on prescribed ABD and shear matrices.
 *
 * The generalized response is evaluated directly from constant section
 * stiffness matrices. If an orientation exists, element strains are rotated
 * into the projected section basis and the response is rotated back into the
 * geometric shell basis required by the element formulation.
 *
 * Physical stress output is an equivalent reconstruction because a prescribed
 * ABD section does not contain a layerwise material model.
 *
 * @see ABDShellSection
 *
 * @author Finn Eggers
 * @date 22.07.2026
 */

#include "section_shell_abd.h"

#include "../core/logging.h"

#include <Eigen/Cholesky>
#include <Eigen/LU>

#include <cmath>
#include <utility>

namespace fem {

/**
 * Constructs and validates a prescribed generalized shell section.
 *
 * The membrane-bending and transverse-shear blocks must contain finite,
 * symmetric and positive-definite coefficients. These requirements ensure a
 * conservative elastic energy and permit the equivalent physical-stress
 * reconstruction used for output. Common thickness and orientation validation
 * is performed by `ShellSection`.
 *
 * @param material Optional material association retained as model metadata.
 * @param region Element region receiving the section.
 * @param thickness Positive physical shell thickness.
 * @param abd Membrane-bending generalized stiffness in the section basis.
 * @param shear Transverse-shear generalized stiffness in the section basis.
 * @param orientation Optional coordinate system defining the section basis.
 * @param csys_axis Zero-based coordinate-system axis projected into the shell plane.
 */
ABDShellSection::ABDShellSection(
    material::Material::Ptr    material,
    model::ElementRegion::Ptr  region,
    Precision                  thickness,
    const Mat6&                abd,
    const Mat2&                shear,
    cos::CoordinateSystem::Ptr orientation,
    Index                      csys_axis
)
    : ShellSection(
          std::move(material),
          std::move(region),
          thickness,
          std::move(orientation),
          csys_axis
      ),
      abd_  (abd),
      shear_(shear) {
    // Finite coefficients are required before symmetry and positive
    // definiteness can be evaluated meaningfully.
    logging::error(abd_.allFinite(),
        "ABDShellSection: ABD matrix must contain only finite values");
    logging::error(shear_.allFinite(),
        "ABDShellSection: shear matrix must contain only finite values");

    // Conservative elastic generalized stiffness matrices are symmetric.
    logging::error(abd_.isApprox(abd_.transpose()),
        "ABDShellSection: ABD matrix must be symmetric");
    logging::error(shear_.isApprox(shear_.transpose()),
        "ABDShellSection: shear matrix must be symmetric");

    // Cholesky factorization checks positive definiteness of both independent
    // generalized stiffness blocks.
    const Eigen::LLT<Mat6> abd_factorization(abd_);
    const Eigen::LLT<Mat2> shear_factorization(shear_);

    logging::error(abd_factorization.info() == Eigen::Success,
        "ABDShellSection: ABD matrix must be positive definite");
    logging::error(shear_factorization.info() == Eigen::Success,
        "ABDShellSection: shear matrix must be positive definite");
}

/**
 * Evaluates prescribed generalized resultants and tangent in the geometric
 * shell basis.
 *
 * Generalized strains are rotated into the projected section basis when an
 * orientation exists. The constant ABD and shear blocks then produce membrane
 * forces, bending moments and transverse shear forces. Resultants and the
 * consistent tangent are transformed back to the geometric shell basis used by
 * element assembly.
 *
 * The formulation is linear elastic and contains no material-point history.
 * Consequently, the state rows and their stride do not change the prescribed
 * response.
 *
 * @param position_reference Physical reference position of the shell point.
 * @param shell_basis_global Geometric shell basis in global coordinates.
 * @param strain_shell Generalized strain in the geometric shell basis.
 * @param old_material_state Unused input state row.
 * @param new_material_state Unused output state row.
 * @param material_state_stride Unused state-row stride.
 * @param resultants_shell Generalized resultants in the geometric shell basis.
 * @param tangent_shell Constant generalized tangent in the geometric shell basis.
 */
void ABDShellSection::evaluate(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    const Precision*              old_material_state,
    Precision*                    new_material_state,
    Index                         material_state_stride,
    ShellStressResultants&        resultants_shell,
    Mat8&                         tangent_shell
) const {
    // A prescribed linear generalized stiffness has no constitutive history.
    (void) old_material_state;
    (void) new_material_state;
    (void) material_state_stride;

    // Without an orientation the prescribed matrices already act in the
    // geometric shell basis. Otherwise construct the projected section basis.
    const Mat3 section_basis_global = orientation_
        ? stress_basis(position_reference, shell_basis_global)
        : shell_basis_global;

    // Express the section in-plane axes in geometric shell coordinates:
    //
    //     R = Q_shell^T Q_section.
    //
    // This rotation maps generalized shell-basis strain components into the
    // section basis used by the prescribed stiffness matrices.
    const Mat2 section_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * section_basis_global.template block<3, 2>(0, 0);

    const Mat8 strain_shell_to_section =
        ShellGeneralizedStrain::transformation(section_axes_in_shell);

    const ShellGeneralizedStrain strain_section(
        strain_shell_to_section * strain_shell.values()
    );

    // Assemble the complete eight-by-eight generalized tangent in the section
    // basis from the prescribed membrane-bending and shear blocks.
    const Index membrane_row = static_cast<Index>(ShellStressResultants::Component::NXX);
    const Index membrane_col = static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);
    const Index shear_row    = static_cast<Index>(ShellStressResultants::Component::QX);
    const Index shear_col    = static_cast<Index>(ShellGeneralizedStrain::Component::GammaXZ);

    Mat8 tangent_section = Mat8::Zero();
    tangent_section.template block<6, 6>(membrane_row, membrane_col) = abd_;
    tangent_section.template block<2, 2>(shear_row, shear_col)       = shear_;

    const ShellStressResultants resultants_section(
        tangent_section * strain_section.values()
    );

    // No transformation is required when geometric shell and section bases are
    // identical.
    if (!orientation_) {
        resultants_shell = resultants_section;
        tangent_shell    = tangent_section;
        return;
    }

    // Physical resultants use a stress-type component transformation. The
    // transpose expresses the geometric shell axes in section coordinates.
    const Mat2 shell_axes_in_section = section_axes_in_shell.transpose();
    const Mat8 resultants_section_to_shell =
        ShellStressResultants::transformation(shell_axes_in_section);

    // Return resultants and tangent in the geometric shell basis used by the B
    // matrix and all element assembly operations.
    resultants_shell = ShellStressResultants(
        resultants_section_to_shell * resultants_section.values()
    );
    tangent_shell = resultants_section_to_shell
                  * tangent_section
                  * strain_shell_to_section;
}

/**
 * Recovers equivalent physical Cauchy stress from an exact base state followed
 * by one affine perturbation.
 *
 * The prescribed ABD law is linear in generalized strain, so the PK2 base
 * stress and PK2 increment are reconstructed directly. The nonlinear
 * PK2-to-Cauchy push-forward is differentiated analytically at F0.
 */
VolumeStressCauchy ABDShellSection::recover_stress(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_base,
    const ShellGeneralizedStrain& strain_increment,
    const Precision*              old_material_state,
    Index                         material_state_stride,
    Precision                     z,
    const Mat3&                   deformation_gradient_base,
    const Mat3&                   deformation_gradient_increment
) const {
    (void) old_material_state;
    (void) material_state_stride;

    const Precision h = thickness_;

    const Mat3 output_basis_global   = stress_basis(position_reference, shell_basis_global);
    const Mat3 recovery_basis_global = orientation_ ? output_basis_global : shell_basis_global;

    const Mat2 recovery_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * recovery_basis_global.template block<3, 2>(0, 0);

    const ShellGeneralizedStrain base_recovery      = strain_base.transformed(recovery_axes_in_shell);
    const ShellGeneralizedStrain increment_recovery = strain_increment.transformed(recovery_axes_in_shell);

    const Index membrane_row = static_cast<Index>(ShellStressResultants::Component::NXX);
    const Index membrane_col = static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);
    const Index shear_row    = static_cast<Index>(ShellStressResultants::Component::QX);
    const Index shear_col    = static_cast<Index>(ShellGeneralizedStrain::Component::GammaXZ);

    Mat8 tangent_recovery = Mat8::Zero();
    tangent_recovery.template block<6, 6>(membrane_row, membrane_col) = abd_;
    tangent_recovery.template block<2, 2>(shear_row, shear_col)       = shear_;

    const ShellStressResultants resultants_base(tangent_recovery * base_recovery.values());
    const ShellStressResultants resultants_increment(tangent_recovery * increment_recovery.values());

    const auto stress_tensor = [h, z](const ShellStressResultants& resultants) {
        const Vec3 plane_stress = resultants.membrane() / h
            + z * (Precision(12) / (h * h * h)) * resultants.moments();
        const Vec2 shear_stress = resultants.transverse_shear() / h;

        Mat3 tensor = Mat3::Zero();
        tensor(0, 0) = plane_stress(0);
        tensor(1, 1) = plane_stress(1);
        tensor(0, 1) = plane_stress(2);
        tensor(1, 0) = plane_stress(2);
        tensor(0, 2) = shear_stress(0);
        tensor(2, 0) = shear_stress(0);
        tensor(1, 2) = shear_stress(1);
        tensor(2, 1) = shear_stress(1);
        return tensor;
    };

    const Mat3 second_pk_base_global =
        recovery_basis_global * stress_tensor(resultants_base) * recovery_basis_global.transpose();
    const Mat3 second_pk_increment_global =
        recovery_basis_global * stress_tensor(resultants_increment) * recovery_basis_global.transpose();

    const Precision J0 = deformation_gradient_base.determinant();
    logging::error(J0 > Precision(0) && std::isfinite(J0),
        "ABDShellSection: invalid base deformation gradient during stress linearization, J = ", J0);

    const Mat3 sigma_base =
        deformation_gradient_base * second_pk_base_global * deformation_gradient_base.transpose() / J0;

    const Mat3 sigma_increment =
        (deformation_gradient_increment * second_pk_base_global * deformation_gradient_base.transpose()
       + deformation_gradient_base * second_pk_increment_global * deformation_gradient_base.transpose()
       + deformation_gradient_base * second_pk_base_global * deformation_gradient_increment.transpose()) / J0
       - (deformation_gradient_base.inverse() * deformation_gradient_increment).trace() * sigma_base;

    const Mat3 sigma = sigma_base + sigma_increment;
    return VolumeStressCauchy(sigma).transformed(Mat3::Identity(), output_basis_global);
}

} // namespace fem
