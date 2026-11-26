/**
 * @file section_shell_composite.h
 * @brief Declares a layered shell section integrated layer by layer.
 *
 * Every layer carries its own material, thickness, in-plane orientation angle
 * and number of Simpson points. Layers are stacked from the bottom surface
 * `z = -h/2` along the shell normal. Generalized strains are reconstructed at
 * the points of every layer, rotated into the layer axes and the material
 * response is integrated into generalized resultants and a consistent tangent.
 *
 * @see CompositeShellSection
 * @see ShellSection
 *
 * @author Finn Eggers
 * @date 24.11.2026
 */

#pragma once

#include "section_shell.h"

#include <string>
#include <vector>

namespace fem {

/**
 * @brief Shell section with explicit layer-wise through-thickness integration.
 *
 * Layer angles rotate the layer material axes about the shell normal, measured
 * from the first section axis. The section axes follow the orientation exactly
 * as for `IntegratedShellSection`; without an orientation they are the
 * geometric shell axes. Material points are enumerated layer by layer from the
 * bottom surface. All layers must share one thermal expansion law because the
 * shell element applies a constant through-thickness free expansion.
 */
struct CompositeShellSection : ShellSection {
    using Ptr = std::shared_ptr<CompositeShellSection>;

    struct Layer {
        fem::material::Material::Ptr material  = nullptr;
        Precision                    thickness = Precision(0);
        Precision                    angle     = Precision(0); // degrees
        Index                        points    = 3;
        std::string                  name;
    };

    // Layers ordered from the bottom to the top shell surface.
    std::vector<Layer> layers_;

    // Construct the section from its layers. The shell thickness is the sum of
    // all layer thicknesses. Every layer needs an elastic material, a positive
    // thickness and an odd number of at least three Simpson points.
    CompositeShellSection(
        std::vector<Layer>         layers,
        model::ElementRegion::Ptr  region,
        cos::CoordinateSystem::Ptr orientation,
        Index                      csys_axis = 0
    );

    // Integrate all layers into membrane forces, moments, corrected transverse
    // shear forces and the consistent generalized tangent. The state pointers
    // identify the first old/new MP rows of the current shell IP.
    void evaluate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        const Precision*              old_material_state,
        Precision*                    new_material_state,
        Index                         material_state_stride,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const override;

    // Recover physical Cauchy stress in the layer containing z from an exact
    // base state followed by one affine perturbation.
    [[nodiscard]] VolumeStressCauchy recover_stress(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_base,
        const ShellGeneralizedStrain& strain_increment,
        const Precision*              old_material_state,
        Index                         material_state_stride,
        Precision                     z,
        const Mat3&                   deformation_gradient_base,
        const Mat3&                   deformation_gradient_increment
    ) const override;

    // Instance copy with remapped region and rigidly placed orientation
    Section::Ptr instantiate(model::ElementRegion::Ptr region,
                             const Mat3&               rotation,
                             const Vec3&               translation) const override;

    // Layer-wise material point lookup, state width and thickness-averaged density
    [[nodiscard]] material::Material::Ptr material_at(Index mp) const override;
    [[nodiscard]] Index                   material_state_size() const override;
    [[nodiscard]] Precision               density() const override;

    // Sum of the Simpson points of all layers
    [[nodiscard]] Index num_mp_per_ip() const override { return num_mp_; }

    void info() override;

private:
    Index num_mp_ = 0;
};

} // namespace fem
