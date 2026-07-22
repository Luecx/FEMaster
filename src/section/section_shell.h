/**
 * @file section_shell.h
 * @brief Declares shell section properties and material matrices.
 *
 * @see src/section/section_shell.cpp
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "section.h"

#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"
#include "../material/strain/shell_generalized_strain.h"
#include "../material/stress/shell_stress_resultants.h"
#include "../material/stress/volume_stress_cauchy.h"

namespace fem {
/**
 * @struct ShellSection
 * @brief Section definition for shell elements.
 */
struct ShellSection : Section {
    using Ptr = std::shared_ptr<ShellSection>; ///< Shared pointer alias for shell sections.

    // Shell thickness and optional material orientation. `csys_axis_` stores
    // the zero-based coordinate-system axis projected into the shell plane.
    Precision                  thickness_   = 1.0;
    cos::CoordinateSystem::Ptr orientation_ = nullptr;
    Index                      csys_axis_   = 0;

    /**
     * @brief Evaluates a shell response in the supplied orthonormal shell basis.
     *
     * The generalized strain is expressed in `shell_basis_global`. Resultants
     * and tangent are returned in the same basis. An optional section
     * orientation is projected into the plane normal to the third basis axis.
     */
    void evaluate(const Vec3&                   position_reference,
                  const Mat3&                   shell_basis_global,
                  const ShellGeneralizedStrain& strain,
                  bool                          use_green_lagrange,
                  ShellStressResultants&        resultants,
                  Mat8&                         tangent) const;

    // Evaluates generalized shell resultants for output. Resultants are
    // returned in the section material basis when an orientation is defined;
    // otherwise a stable global-axis projection is used.
    [[nodiscard]] ShellStressResultants compute_resultants(const Vec3&                   position_reference,
                                                           const Mat3&                   shell_basis_global,
                                                           const ShellGeneralizedStrain& strain,
                                                           bool                          use_green_lagrange) const;

    // Evaluates a reconstructed through-thickness Cauchy stress for output.
    // Without an orientation the stress is returned in global components; with
    // an orientation it is returned in the projected section material basis.
    [[nodiscard]] virtual VolumeStressCauchy compute_stress(const Vec3&                   position_reference,
                                                            const Mat3&                   shell_basis_global,
                                                            const ShellGeneralizedStrain& strain,
                                                            Precision                     z,
                                                            bool                          use_green_lagrange,
                                                            Mat3                          deformation_gradient = Mat3::Identity()) const;

    virtual Index num_mp_per_ip() const { return 1; }

    /**
     * @brief Outputs shell section details through the logger.
     */
    void info() override;

    /**
     * @brief Builds a compact one-line shell section summary.
     *
     * @return std::string Material, region, density, orientation and thickness.
     */
    std::string str() const override;

protected:
    ShellSection(material::Material::Ptr    material,
                 model::ElementRegion::Ptr  region,
                 Precision                  thickness,
                 cos::CoordinateSystem::Ptr orientation,
                 Index                      csys_axis = 0);

    // Output and material basis construction. The returned basis has the
    // selected or fallback direction as first in-plane axis and the supplied
    // shell normal as third axis.
    [[nodiscard]] Mat3 output_basis_global(const Vec3& position_reference,
                                           const Mat3& shell_basis_global) const;

    [[nodiscard]] Mat2 output_to_shell_rotation(const Vec3& position_reference,
                                                const Mat3& shell_basis_global) const;

    virtual void evaluate_material(const ShellGeneralizedStrain& strain,
                                   bool                          use_green_lagrange,
                                   ShellStressResultants&        resultants,
                                   Mat8&                         tangent) const = 0;
};

} // namespace fem
