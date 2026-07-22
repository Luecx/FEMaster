/**
 * @file section_shell.h
 * @brief Declares the abstract base class shared by all shell section formulations.
 *
 * A shell element supplies generalized strains in its pointwise geometric shell
 * basis. Concrete section formulations convert these strains into generalized
 * membrane forces, bending moments, transverse shear forces and a consistent
 * section tangent.
 *
 * `ShellSection` does not implement a constitutive law itself. It owns only the
 * data and output conventions shared by every shell section: thickness,
 * optional orientation, selected coordinate-system axis and the bases used for
 * physical stress and generalized stress-resultant output.
 *
 * The element-facing `evaluate()` contract is identical for every section:
 * input strains, output resultants and the tangent are all expressed in the
 * supplied geometric shell basis. Concrete sections are responsible for any
 * temporary transformation into their material basis.
 *
 * @see ABDShellSection
 * @see IntegratedShellSection
 *
 * @author Finn Eggers
 * @date 22.07.2026
 */

#pragma once

#include "section.h"

#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"
#include "../material/strain/shell_generalized_strain.h"
#include "../material/stress/shell_stress_resultants.h"
#include "../material/stress/volume_stress_cauchy.h"

#include <memory>
#include <string>

namespace fem {

/**
 * @brief Abstract common base for shell section formulations.
 *
 * The class stores shared section properties and implements only behavior that
 * is independent of the constitutive formulation. It intentionally contains no
 * default generalized response and no intermediate virtual `evaluate_material`
 * hook. Every concrete section implements `evaluate()` and physical stress
 * recovery directly.
 *
 * An optional coordinate system defines both material and local output axes.
 * `csys_axis_` selects the zero-based coordinate-system axis projected into the
 * shell plane. A nearly vanishing projection is a hard model error; an explicit
 * orientation is never replaced silently by an element-local direction.
 */
struct ShellSection : Section {
    using Ptr = std::shared_ptr<ShellSection>;

    // Physical shell thickness used by constitutive integration, stress
    // recovery, mass integration and through-thickness coordinate mappings.
    Precision thickness_ = Precision(1);

    // Optional spatial coordinate system defining material and output axes.
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Zero-based axis selected by the external one-based CSYSAXIS convention.
    Index csys_axis_ = 0;

    ~ShellSection() override = default;

    /**
     * @brief Evaluates generalized resultants and the consistent section tangent.
     *
     * @param position_reference Physical reference position of the evaluation point.
     * @param shell_basis_global Orthonormal geometric shell basis in global coordinates.
     * @param strain_shell Generalized strain in the geometric shell basis.
     * @param use_green_lagrange Select finite-strain material evaluation.
     * @param resultants_shell Resultants returned in the geometric shell basis.
     * @param tangent_shell Tangent returned in the geometric shell basis.
     */
    virtual void evaluate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        bool                          use_green_lagrange,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const = 0;

    /**
     * @brief Evaluates generalized resultants in the configured output basis.
     *
     * The concrete section response is first evaluated in the geometric shell
     * basis through `evaluate()`. The common base implementation then rotates
     * the resultants into `stress_resultant_basis()`.
     *
     * @param position_reference Physical reference position of the evaluation point.
     * @param shell_basis_global Orthonormal geometric shell basis in global coordinates.
     * @param strain_shell Generalized strain in the geometric shell basis.
     * @param use_green_lagrange Select finite-strain material evaluation.
     * @return Resultants in the configured stress-resultant output basis.
     */
    [[nodiscard]] ShellStressResultants evaluate_output_resultants(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        bool                          use_green_lagrange
    ) const;

    /**
     * @brief Evaluates physical Cauchy stress at one thickness coordinate.
     *
     * Without an explicit orientation, stress components are global. With an
     * orientation, they are expressed in the projected section basis.
     *
     * @param position_reference Physical reference position of the evaluation point.
     * @param shell_basis_global Orthonormal geometric shell basis in global coordinates.
     * @param strain_shell Generalized strain in the geometric shell basis.
     * @param z Physical thickness coordinate measured from the midsurface.
     * @param use_green_lagrange Select finite-strain PK2-to-Cauchy recovery.
     * @param deformation_gradient Three-dimensional deformation gradient.
     * @return Cauchy stress in the configured physical-stress basis.
     */
    [[nodiscard]] virtual VolumeStressCauchy evaluate_output_stress(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        Precision                     z,
        bool                          use_green_lagrange,
        const Mat3&                   deformation_gradient = Mat3::Identity()
    ) const = 0;

    /**
     * @brief Returns the number of material points per shell integration point.
     *
     * @return Number of through-thickness material points.
     */
    [[nodiscard]] virtual Index num_mp_per_ip() const = 0;

    // Output the common section properties through the project logger.
    void info() override;

    // Build a compact one-line representation of the common section properties.
    [[nodiscard]] std::string str() const override;

protected:
    /**
     * @brief Initializes common shell-section data.
     *
     * @param material Material assigned to the section.
     * @param region Element region receiving the section.
     * @param thickness Positive physical shell thickness.
     * @param orientation Optional coordinate system defining section directions.
     * @param csys_axis Zero-based selected coordinate-system axis.
     */
    ShellSection(
        material::Material::Ptr    material,
        model::ElementRegion::Ptr  region,
        Precision                  thickness,
        cos::CoordinateSystem::Ptr orientation,
        Index                      csys_axis = 0
    );

    /**
     * @brief Builds the basis used for physical shell stress output.
     *
     * Without an orientation, the global Cartesian basis is returned. With an
     * orientation, the selected coordinate-system axis is projected into the
     * shell plane and completed with the shell normal to a right-handed basis.
     *
     * @param position_reference Physical reference position of the evaluation point.
     * @param shell_basis_global Orthonormal shell basis whose third axis is the normal.
     * @return Stress basis expressed in global coordinates.
     */
    [[nodiscard]] Mat3 stress_basis(
        const Vec3& position_reference,
        const Mat3& shell_basis_global
    ) const;

    /**
     * @brief Builds the basis used for generalized stress-resultant output.
     *
     * With an orientation, this equals `stress_basis()`. Without one, global X
     * is projected into the shell plane and global Y is used only if required.
     *
     * @param position_reference Physical reference position of the evaluation point.
     * @param shell_basis_global Orthonormal geometric shell basis.
     * @return Stress-resultant basis expressed in global coordinates.
     */
    [[nodiscard]] Mat3 stress_resultant_basis(
        const Vec3& position_reference,
        const Mat3& shell_basis_global
    ) const;
};

} // namespace fem
