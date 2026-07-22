/**
 * @file section_shell_abd.h
 * @brief Declares the shell section based on prescribed ABD and shear matrices.
 *
 * @see src/section/section_shell_abd.cpp
 * @author Finn Eggers
 * @date 18.07.2026
 */

#pragma once

#include "section_shell.h"

namespace fem {
struct ABDShellSection : ShellSection {
    using Ptr = std::shared_ptr<ABDShellSection>;

    Mat6 abd_   = Mat6::Zero();
    Mat2 shear_ = Mat2::Zero();

    ABDShellSection(material::Material::Ptr    material,
                    model::ElementRegion::Ptr  region,
                    Precision                  thickness,
                    const Mat6&                abd,
                    const Mat2&                shear,
                    cos::CoordinateSystem::Ptr orientation,
                    Index                      csys_axis = 0);

    Index num_mp_per_ip() const override { return 0; }

    [[nodiscard]] VolumeStressCauchy compute_stress(const Vec3&                   position_reference,
                                                    const Mat3&                   shell_basis_global,
                                                    const ShellGeneralizedStrain& strain,
                                                    Precision                     z,
                                                    bool                          use_green_lagrange,
                                                    Mat3                          deformation_gradient = Mat3::Identity()) const override;

    void info() override;
    std::string str() const override;

protected:
    void evaluate_material(const ShellGeneralizedStrain& strain,
                           bool                          use_green_lagrange,
                           ShellStressResultants&        resultants,
                           Mat8&                         tangent) const override;
};
} // namespace fem
