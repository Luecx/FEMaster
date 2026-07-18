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

namespace fem {
/**
 * @struct ShellSection
 * @brief Section definition for shell elements.
 */
struct ShellSection : Section {
    using Ptr = std::shared_ptr<ShellSection>; ///< Shared pointer alias for shell sections.

    Precision                  thickness_   = 1.0;     ///< Shell thickness.
    cos::CoordinateSystem::Ptr orientation_ = nullptr; ///< Optional section/material orientation.

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
                 cos::CoordinateSystem::Ptr orientation);

    virtual void evaluate_material(const ShellGeneralizedStrain& strain,
                                   bool                          use_green_lagrange,
                                   ShellStressResultants&        resultants,
                                   Mat8&                         tangent) const = 0;
};

struct IntegratedShellSection : ShellSection {
    using Ptr = std::shared_ptr<IntegratedShellSection>;

    IntegratedShellSection(material::Material::Ptr    material,
                           model::ElementRegion::Ptr  region,
                           Precision                  thickness,
                           cos::CoordinateSystem::Ptr orientation);

    Index num_mp_per_ip() const override { return 5; }

protected:
    void evaluate_material(const ShellGeneralizedStrain& strain,
                           bool                          use_green_lagrange,
                           ShellStressResultants&        resultants,
                           Mat8&                         tangent) const override;
};

struct ABDShellSection : ShellSection {
    using Ptr = std::shared_ptr<ABDShellSection>;

    Mat6 abd_   = Mat6::Zero();
    Mat2 shear_ = Mat2::Zero();

    ABDShellSection(material::Material::Ptr    material,
                    model::ElementRegion::Ptr  region,
                    Precision                  thickness,
                    const Mat6&                abd,
                    const Mat2&                shear,
                    cos::CoordinateSystem::Ptr orientation);

    Index num_mp_per_ip() const override { return 0; }

    void info() override;
    std::string str() const override;

protected:
    void evaluate_material(const ShellGeneralizedStrain& strain,
                           bool                          use_green_lagrange,
                           ShellStressResultants&        resultants,
                           Mat8&                         tangent) const override;
};
} // namespace fem
