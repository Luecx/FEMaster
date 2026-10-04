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
 * data shared by every shell section: thickness, optional orientation and the
 * selected coordinate-system axis used for constitutive and physical-stress recovery.
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
 * recovery directly through recover_stress().
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
    // outputs use the supplied geometric shell basis. old_material_state and
    // new_material_state address the first through-thickness input/output rows
    // at the current shell IP; material_state_stride is the scalar distance to
    // each following row.
    // Concrete sections perform all material-basis transformations internally.
    virtual void evaluate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        const Precision*              old_material_state,
        Precision*                    new_material_state,
        Index                         material_state_stride,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const = 0;

    // Recover physical Cauchy stress from an exact base state followed by
    // one affine perturbation. strain_increment is the mechanical generalized
    // strain increment from the base state and deformation_gradient_increment
    // is the matching first variation of F. Constitutive history remains read-only.
    [[nodiscard]] virtual VolumeStressCauchy recover_stress(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_base,
        const ShellGeneralizedStrain& strain_increment,
        const Precision*              old_material_state,
        Index                         material_state_stride,
        Precision                     z,
        const Mat3&                   deformation_gradient_base,
        const Mat3&                   deformation_gradient_increment
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

};

} // namespace fem
