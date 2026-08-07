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

    // Constitutive evaluation in the global solid basis
    void evaluate(const Vec3&                   position_reference,
                  const Mat3&                   additional_rotation,
                  const VolumeStrainLinearized& strain_global,
                  Precision*                    state,
                  VolumeStressCauchy&           stress_global,
                  Mat6&                         tangent_global) const;

    void evaluate(const Vec3&                      position_reference,
                  const Mat3&                      additional_rotation,
                  const VolumeStrainGreenLagrange& strain_global,
                  Precision*                       state,
                  VolumeStressPK2&                 stress_global,
                  Mat6&                            tangent_global) const;

    // Tangent derivatives with respect to the three additional material rotations
    std::array<Mat6, 3> tangent_rotation_derivatives(
        const Vec3&                position_reference,
        const Mat3&                additional_rotation,
        const std::array<Mat3, 3>& additional_rotation_derivatives,
        Precision*                 state
    ) const;

    // Section diagnostics
    void info() override;
    std::string str() const override;

private:
    // Material basis at one reference position
    Mat3 section_orientation_basis(const Vec3& position_reference) const;
};

} // namespace fem
