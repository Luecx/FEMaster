/**
 * @file section_solid.h
 * @brief Declares the constitutive section used by three-dimensional solids.
 *
 * @author Finn Eggers
 * @date 16.08.2026
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

    void update(const Vec3& position_reference,
                const Mat3& additional_rotation,
                const VolumeStrainLinearized& strain_global,
                const Precision* state_old,
                Precision* state_new,
                VolumeStressCauchy& stress_global,
                Mat6& tangent_global) const;

    void update(const Vec3& position_reference,
                const Mat3& additional_rotation,
                const VolumeStrainGreenLagrange& strain_global,
                const Precision* state_old,
                Precision* state_new,
                VolumeStressPK2& stress_global,
                Mat6& tangent_global) const;

    void recover(const Vec3& position_reference,
                 const Mat3& additional_rotation,
                 const VolumeStrainLinearized& strain_global,
                 const Precision* state,
                 VolumeStressCauchy& stress_global) const;

    void recover(const Vec3& position_reference,
                 const Mat3& additional_rotation,
                 const VolumeStrainGreenLagrange& strain_global,
                 const Precision* state,
                 VolumeStressPK2& stress_global) const;

    // State-neutral tangent of the elastic backbone at the undeformed state.
    // Auxiliary formulations such as hourglass control must use this query
    // instead of initiating a plastic constitutive update.
    Mat6 elastic_tangent_reference(const Vec3& position_reference,
                                   const Mat3& additional_rotation) const;

    std::array<Mat6, 3> tangent_rotation_derivatives(
        const Vec3& position_reference,
        const Mat3& additional_rotation,
        const std::array<Mat3, 3>& additional_rotation_derivatives
    ) const;

    void info() override;
    std::string str() const override;

private:
    Mat3 section_orientation_basis(const Vec3& position_reference) const;
};

} // namespace fem
