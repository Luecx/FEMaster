/**
 * @file section_shell_abd.h
 * @brief Declares a shell section defined by prescribed ABD and shear matrices.
 *
 * The section contains no through-thickness material-point state. Its generalized
 * response follows directly from constant section stiffness matrices, while
 * physical stress output is reconstructed as one equivalent homogeneous layer.
 *
 * @see ABDShellSection
 * @see ShellSection
 *
 * @author Finn Eggers
 * @date 22.07.2026
 */

#pragma once

#include "section_shell.h"

namespace fem {

/**
 * @brief Shell section with prescribed generalized stiffness matrices.
 *
 * `abd_` maps membrane strains and curvatures to membrane forces and bending
 * moments. `shear_` maps transverse engineering shear strains to transverse
 * shear forces. Both matrices are interpreted in the user-defined section basis
 * when an orientation exists and in the geometric shell basis otherwise.
 */
struct ABDShellSection : ShellSection {
    using Ptr = std::shared_ptr<ABDShellSection>;

    // Membrane-bending stiffness matrix in the section basis.
    Mat6 abd_ = Mat6::Zero();

    // Transverse-shear stiffness matrix in the section basis.
    Mat2 shear_ = Mat2::Zero();

    /**
     * @brief Constructs and validates a prescribed generalized shell section.
     *
     * @param material Optional material association retained as model metadata.
     * @param region Element region receiving the section.
     * @param thickness Positive physical shell thickness.
     * @param abd Symmetric positive-definite membrane-bending stiffness matrix.
     * @param shear Symmetric positive-definite transverse-shear stiffness matrix.
     * @param orientation Optional coordinate system defining the section basis.
     * @param csys_axis Zero-based coordinate-system axis projected into the shell plane.
     */
    ABDShellSection(
        material::Material::Ptr    material,
        model::ElementRegion::Ptr  region,
        Precision                  thickness,
        const Mat6&                abd,
        const Mat2&                shear,
        cos::CoordinateSystem::Ptr orientation,
        Index                      csys_axis = 0
    );

    /**
     * @brief Evaluates prescribed generalized resultants in the shell basis.
     *
     * @copydetails ShellSection::evaluate
     */
    void evaluate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        bool                          use_green_lagrange,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const override;

    /**
     * @brief Reconstructs an equivalent through-thickness Cauchy stress.
     *
     * @copydetails ShellSection::evaluate_output_stress
     */
    [[nodiscard]] VolumeStressCauchy evaluate_output_stress(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        Precision                     z,
        bool                          use_green_lagrange,
        const Mat3&                   deformation_gradient
    ) const override;

    // A prescribed ABD section contains no actual material points.
    [[nodiscard]] Index num_mp_per_ip() const override { return 0; }
};

} // namespace fem
