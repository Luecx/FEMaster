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

    // Construct a section whose material response is sampled at five composite-
    // Simpson points through the physical thickness. orientation and csys_axis
    // define the local material basis used consistently at all five points.
    IntegratedShellSection(
        material::Material::Ptr    material,
        model::ElementRegion::Ptr  region,
        Precision                  thickness,
        cos::CoordinateSystem::Ptr orientation,
        Index                      csys_axis = 0
    );

    // Reconstruct epsilon(z) = epsilon_0 + z kappa at five physical material
    // points and integrate stress into membrane forces, moments and corrected
    // transverse shear forces. material_state identifies the first MP row and
    // material_state_stride advances through its four following rows. Each row
    // is passed directly to the selected linearized or Green-Lagrange material
    // evaluation, and the consistent tangent is integrated by the same rule.
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

    // Recover physical Cauchy stress at arbitrary z using the closest stored
    // Simpson material point as its history state. Strain is reconstructed in
    // the section material basis. PK2 output from Green-Lagrange evaluation is
    // pushed forward with deformation_gradient; linearized output is already
    // Cauchy stress. The final components follow the common output-basis convention.
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

    // Five composite-Simpson material points are enumerated for every in-plane
    // shell integration point and stored as consecutive state rows.
    [[nodiscard]] Index num_mp_per_ip() const override { return 5; }
};

} // namespace fem
