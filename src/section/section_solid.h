/**
 * @file section_solid.h
 * @brief Declares the constitutive section used by three-dimensional solid elements.
 *
 * `SolidSection` transforms global solid strains into the optional material
 * coordinate system, evaluates the assigned elasticity model and transforms
 * stresses and optional consistent tangents back into the global basis. The
 * element owns material-point addressing and supplies separate input and output
 * state rows.
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
 * history. Every constitutive call receives the old and new state rows selected
 * by the element. Linearized response returns Cauchy stress, while finite-strain
 * response uses Green-Lagrange strain and PK2 stress in the reference material
 * basis before transformation back to global coordinates.
 *
 * Tangent output is optional. A null tangent pointer is forwarded to the material
 * model so stress-only nonlinear residual evaluations can avoid constitutive and
 * coordinate-transformation work that is required only for stiffness assembly.
 */
struct SolidSection : Section {
    using Ptr = std::shared_ptr<SolidSection>;

    // Optional material orientation in the reference configuration
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Evaluate infinitesimal solid response. Global engineering strain is
    // transformed into the material basis and Cauchy stress is transformed back.
    // The globally transformed tangent is produced only for a non-null pointer.
    void evaluate(const Vec3&                   position_reference,
                  const Mat3&                   additional_rotation,
                  const VolumeStrainLinearized& strain_global,
                  const Precision*              old_state,
                  Precision*                    new_state,
                  VolumeStressCauchy&           stress_global,
                  Mat6*                         tangent_global = nullptr) const;

    // Evaluate Total-Lagrangian solid response with the same basis and state
    // ownership convention. Green-Lagrange strain and PK2 stress remain in the
    // reference configuration; tangent output is optional.
    void evaluate(const Vec3&                      position_reference,
                  const Mat3&                      additional_rotation,
                  const VolumeStrainGreenLagrange& strain_global,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress_global,
                  Mat6*                            tangent_global = nullptr) const;

    // Differentiate the globally transformed linear material tangent with
    // respect to three supplied additional-rotation directions. This operation
    // explicitly requires the material tangent at zero strain.
    std::array<Mat6, 3> tangent_rotation_derivatives(
        const Vec3&                position_reference,
        const Mat3&                additional_rotation,
        const std::array<Mat3, 3>& additional_rotation_derivatives,
        const Precision*           old_state,
        Precision*                 new_state
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
