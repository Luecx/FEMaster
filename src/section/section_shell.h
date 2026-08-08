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

    // Lifetime through the polymorphic section interface
    ~ShellSection() override = default;

    // Evaluate generalized membrane forces, bending moments, transverse shear
    // forces and their consistent eight-by-eight tangent. Input strain and both
    // outputs use the supplied geometric shell basis. material_state addresses
    // the first through-thickness material point at the current shell IP;
    // material_state_stride is the scalar distance to each following row.
    // Concrete sections perform all material-basis transformations internally.
    virtual void evaluate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        Precision*                    material_state,
        Index                         material_state_stride,
        bool                          use_green_lagrange,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const = 0;

    // Recover generalized resultants for output. The concrete section is first
    // evaluated with the same state pointer, stride and strain-measure contract
    // as evaluate(). The base then rotates physical membrane, moment and shear
    // components from the geometric shell basis into stress_resultant_basis().
    // Any in-place constitutive state update occurs during the concrete call.
    [[nodiscard]] ShellStressResultants evaluate_output_resultants(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        Precision*                    material_state,
        Index                         material_state_stride,
        bool                          use_green_lagrange
    ) const;

    // Recover physical Cauchy stress at thickness coordinate z measured from
    // the midsurface. Linearized sections return Cauchy stress directly;
    // finite-strain sections use deformation_gradient to push their PK2 material
    // response forward. Components are global without an orientation and
    // section-local with one. material_state identifies the first state row at
    // the parent shell IP and the concrete formulation chooses the relevant MP.
    [[nodiscard]] virtual VolumeStressCauchy evaluate_output_stress(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        Precision*                    material_state,
        Index                         material_state_stride,
        Precision                     z,
        bool                          use_green_lagrange,
        const Mat3&                   deformation_gradient = Mat3::Identity()
    ) const = 0;

    // Return the fixed number of constitutive material points stored for every
    // in-plane shell integration point. Element MP enumeration uses this count
    // to allocate a contiguous state-row block before any constitutive call.
    [[nodiscard]] virtual Index num_mp_per_ip() const = 0;

    // Output material, region, orientation, selected axis and thickness through
    // the project logger using the common shell-section representation.
    void info() override;

    // Build a stable one-line summary of the same common section properties for
    // model diagnostics and stream output.
    [[nodiscard]] std::string str() const override;

protected:
    // Initialize and validate common section data. Thickness must be positive;
    // csys_axis is the internal zero-based axis projected into the shell plane.
    // Construction is restricted to concrete shell-section formulations.
    ShellSection(
        material::Material::Ptr    material,
        model::ElementRegion::Ptr  region,
        Precision                  thickness,
        cos::CoordinateSystem::Ptr orientation,
        Index                      csys_axis = 0
    );

    // Build the physical-stress output basis in global coordinates. Without an
    // orientation this is the global Cartesian basis. Otherwise the selected
    // coordinate-system axis is projected into the shell tangent plane and
    // completed with the supplied shell normal to a right-handed basis.
    [[nodiscard]] Mat3 stress_basis(
        const Vec3& position_reference,
        const Mat3& shell_basis_global
    ) const;

    // Build the generalized-resultant output basis in global coordinates. An
    // explicit orientation uses the physical-stress basis. Without one, global
    // X, or global Y as a geometric fallback, defines a deterministic tangent
    // axis for component-wise averaging of neighboring shell resultants.
    [[nodiscard]] Mat3 stress_resultant_basis(
        const Vec3& position_reference,
        const Mat3& shell_basis_global
    ) const;
};

} // namespace fem
