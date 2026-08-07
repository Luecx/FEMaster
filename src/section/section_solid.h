/**
 * @file section_solid.h
 * @brief Declares the constitutive section used by three-dimensional solid elements.
 *
 * `SolidSection` transforms global solid strains into the optional material
 * coordinate system, evaluates the assigned elasticity model and transforms
 * stresses and consistent tangents back into the global basis. The element owns
 * material-point addressing and supplies the active state pointer directly.
 *
 * @see SolidSection
 * @see src/section/section_solid.cpp
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "section.h"

#include "../cos/coordinate_system.h"
#include "../material/strain/volume_strain_green_lagrange.h"
#include "../material/strain/volume_strain_linearized.h"
#include "../material/stress/volume_stress_cauchy.h"
#include "../material/stress/volume_stress_pk2.h"

#include <array>

namespace fem {

/**
 * @brief Constitutive section for three-dimensional solid elements.
 *
 * The section owns the optional material orientation but no material-point
 * history. Every constitutive call receives the state row selected by the
 * element. Linearized response returns Cauchy stress, while finite-strain
 * response uses Green-Lagrange strain and PK2 stress in the reference material
 * basis before transformation back to global coordinates.
 */
struct SolidSection : Section {
    using Ptr = std::shared_ptr<SolidSection>;

    // Optional material orientation in the reference configuration
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Evaluate infinitesimal solid response. Global engineering strain is
    // transformed into the reference material basis composed from the optional
    // section orientation and additional element rotation. The selected state
    // row is passed directly to the material. Cauchy stress and the consistently
    // transformed six-by-six tangent are returned in global coordinates.
    void evaluate(const Vec3&                   position_reference,
                  const Mat3&                   additional_rotation,
                  const VolumeStrainLinearized& strain_global,
                  Precision*                    state,
                  VolumeStressCauchy&           stress_global,
                  Mat6&                         tangent_global) const;

    // Evaluate Total-Lagrangian solid response with the same basis and state
    // ownership convention. Green-Lagrange strain is transformed into material
    // coordinates; second Piola-Kirchhoff stress and dS/dE are transformed back
    // into the global reference basis without changing their stress measure.
    void evaluate(const Vec3&                      position_reference,
                  const Mat3&                      additional_rotation,
                  const VolumeStrainGreenLagrange& strain_global,
                  Precision*                       state,
                  VolumeStressPK2&                 stress_global,
                  Mat6&                            tangent_global) const;

    // Differentiate the globally transformed linear material tangent with
    // respect to three supplied additional-rotation directions. The material
    // tangent is evaluated at zero strain using the active state row and is
    // assumed independent of the rotation parameters; only strain and stress
    // transformation operators are differentiated.
    std::array<Mat6, 3> tangent_rotation_derivatives(
        const Vec3&                position_reference,
        const Mat3&                additional_rotation,
        const std::array<Mat3, 3>& additional_rotation_derivatives,
        Precision*                 state
    ) const;

    // Report material, region and optional orientation through the logger and a
    // compact stable string representation.
    void info() override;
    std::string str() const override;

private:
    // Evaluate the optional spatial coordinate system at one physical reference
    // point. Without an orientation, return the global Cartesian basis.
    Mat3 section_orientation_basis(const Vec3& position_reference) const;
};

} // namespace fem
