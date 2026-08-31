/**
 * @file section_solid.h
 * @brief Declares the constitutive section used by three-dimensional solid elements.
 *
 * `SolidSection` transforms global solid strains into the optional material
 * coordinate system and dispatches either read-only constitutive evaluation or
 * explicit source-to-target history integration.
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

struct SolidSection : Section {
    using Ptr = std::shared_ptr<SolidSection>;

    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Read-only response used by recovery, diagnostics and state-neutral paths.
    void evaluate(const Vec3&                   position_reference,
                  const Mat3&                   additional_rotation,
                  const VolumeStrainLinearized& strain_global,
                  const Precision*              state,
                  VolumeStressCauchy&           stress_global,
                  Mat6&                         tangent_global) const;

    void evaluate(const Vec3&                      position_reference,
                  const Mat3&                      additional_rotation,
                  const VolumeStrainGreenLagrange& strain_global,
                  const Precision*                 state,
                  VolumeStressPK2&                 stress_global,
                  Mat6&                            tangent_global) const;

    // Constitutive integration used by nonlinear assembly. Source history is
    // immutable and the complete updated history is written to target_state.
    void integrate(const Vec3&                   position_reference,
                   const Mat3&                   additional_rotation,
                   const VolumeStrainLinearized& strain_global,
                   const Precision*              state,
                   Precision*                    target_state,
                   VolumeStressCauchy&           stress_global,
                   Mat6&                         tangent_global) const;

    void integrate(const Vec3&                      position_reference,
                   const Mat3&                      additional_rotation,
                   const VolumeStrainGreenLagrange& strain_global,
                   const Precision*                 state,
                   Precision*                       target_state,
                   VolumeStressPK2&                 stress_global,
                   Mat6&                            tangent_global) const;

    std::array<Mat6, 3> tangent_rotation_derivatives(
        const Vec3&                position_reference,
        const Mat3&                additional_rotation,
        const std::array<Mat3, 3>& additional_rotation_derivatives,
        const Precision*           state
    ) const;

    void info() override;
    std::string str() const override;

private:
    Mat3 section_orientation_basis(const Vec3& position_reference) const;
};

} // namespace fem
