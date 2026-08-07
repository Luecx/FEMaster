/**
 * @file section_shell_abd.h
 * @brief Declares a shell section defined by prescribed ABD and shear matrices.
 *
 * The generalized response follows directly from constant section stiffness
 * matrices, while physical stress output is reconstructed as one equivalent
 * homogeneous layer. One material-state slot per shell integration point keeps
 * the common element-to-section state interface uniform; the prescribed ABD
 * formulation itself does not modify that state.
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

    // Construct a prescribed elastic section and validate that both generalized
    // stiffness blocks are finite, symmetric and positive definite. thickness
    // controls equivalent physical-stress recovery; orientation and csys_axis
    // define the basis in which abd and shear are interpreted.
    ABDShellSection(
        material::Material::Ptr    material,
        model::ElementRegion::Ptr  region,
        Precision                  thickness,
        const Mat6&                abd,
        const Mat2&                shear,
        cos::CoordinateSystem::Ptr orientation,
        Index                      csys_axis = 0
    );

    // Apply the constant ABD and transverse-shear operators. Generalized strain
    // is rotated into the prescribed section basis and resultants/tangent are
    // returned in the geometric shell basis. The formulation has no history:
    // material_state and material_state_stride are accepted for interface
    // uniformity but never read or modified, and the strain-measure selector
    // does not alter the linear generalized law.
    void evaluate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        Precision*                    material_state,
        Index                         material_state_stride,
        bool                          use_green_lagrange,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const override;

    // Reconstruct one equivalent homogeneous-layer stress distribution from
    // N/h + 12 z M/h^3 and Q/h. Linearized recovery is already Cauchy stress;
    // finite-strain recovery treats the reconstructed tensor as PK2 stress and
    // pushes it forward with deformation_gradient. Output follows the common
    // global/section-local basis convention and remains state-neutral.
    [[nodiscard]] VolumeStressCauchy evaluate_output_stress(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        Precision*                    material_state,
        Index                         material_state_stride,
        Precision                     z,
        bool                          use_green_lagrange,
        const Mat3&                   deformation_gradient
    ) const override;

    // Reserve one state row per shell IP so element material-point addressing is
    // uniform across section types. The prescribed ABD formulation does not use it.
    [[nodiscard]] Index num_mp_per_ip() const override { return 1; }
};

} // namespace fem
