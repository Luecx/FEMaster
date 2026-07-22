/**
 * @file section_shell_integrated.h
 * @brief Declares a shell section integrated through the physical thickness.
 *
 * Generalized shell strains are reconstructed at five Simpson points through
 * the thickness. The assigned material model is evaluated at every point and
 * the resulting stresses and tangents are integrated into generalized shell
 * resultants and a consistent section tangent.
 *
 * @see IntegratedShellSection
 * @see ShellSection
 *
 * @author Finn Eggers
 * @date 22.07.2026
 */

#pragma once

#include "section_shell.h"

namespace fem {

/**
 * @brief Shell section with explicit through-thickness material integration.
 *
 * Without an orientation, the geometric shell basis is also the material basis.
 * With an orientation, generalized strains are rotated into the projected
 * section basis before material evaluation and the integrated response is
 * rotated back into the geometric shell basis for assembly.
 */
struct IntegratedShellSection : ShellSection {
    using Ptr = std::shared_ptr<IntegratedShellSection>;

    /**
     * @brief Constructs a through-thickness integrated shell section.
     *
     * @param material Material providing shell-point elasticity evaluation.
     * @param region Element region receiving the section.
     * @param thickness Positive physical shell thickness.
     * @param orientation Optional coordinate system defining the section basis.
     * @param csys_axis Zero-based coordinate-system axis projected into the shell plane.
     */
    IntegratedShellSection(
        material::Material::Ptr    material,
        model::ElementRegion::Ptr  region,
        Precision                  thickness,
        cos::CoordinateSystem::Ptr orientation,
        Index                      csys_axis = 0
    );

    /**
     * @brief Integrates generalized resultants and tangent through the thickness.
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
     * @brief Evaluates material stress at one physical thickness coordinate.
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

    // Five composite Simpson points are used for every shell integration point.
    [[nodiscard]] Index num_mp_per_ip() const override { return 5; }
};

} // namespace fem
