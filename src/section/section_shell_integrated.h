/**
 * @file section_shell_integrated.h
 * @brief Declares the through-thickness integrated shell section.
 *
 * @see src/section/section_shell_integrated.cpp
 * @author Finn Eggers
 * @date 18.07.2026
 */

#pragma once

#include "section_shell.h"

namespace fem {
struct IntegratedShellSection : ShellSection {
    using Ptr = std::shared_ptr<IntegratedShellSection>;

    IntegratedShellSection(material::Material::Ptr    material,
                           model::ElementRegion::Ptr  region,
                           Precision                  thickness,
                           cos::CoordinateSystem::Ptr orientation,
                           Index                      csys_axis = 0);

    Index num_mp_per_ip() const override { return 5; }

    [[nodiscard]] VolumeStressCauchy compute_stress(const Vec3&                   position_reference,
                                                    const Mat3&                   shell_basis_global,
                                                    const ShellGeneralizedStrain& strain,
                                                    Precision                     z,
                                                    bool                          use_green_lagrange,
                                                    Mat3                          deformation_gradient = Mat3::Identity()) const override;

protected:
    void evaluate_material(const ShellGeneralizedStrain& strain,
                           bool                          use_green_lagrange,
                           ShellStressResultants&        resultants,
                           Mat8&                         tangent) const override;
};
} // namespace fem
