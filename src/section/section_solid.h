/**
 * @file section_solid.h
 * @brief Declares solid section properties.
 *
 * @see src/section/section_solid.cpp
 * @author Finn Eggers
 * @date 28.04.2026
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
 * @struct SolidSection
 * @brief Section definition for solid elements.
 */
struct SolidSection : Section {
    using Ptr = std::shared_ptr<SolidSection>; ///< Shared pointer alias for solid sections.

    cos::CoordinateSystem::Ptr orientation_ = nullptr; ///< Optional material orientation.

    /**
     * @brief Evaluates a linearized solid response in the global basis.
     */
    void evaluate(const Vec3&                   position_reference,
                  const Mat3&                   additional_rotation,
                  const VolumeStrainLinearized& strain_global,
                  VolumeStressCauchy&            stress_global,
                  Mat6&                          tangent_global) const;

    /**
     * @brief Evaluates a Green-Lagrange solid response in the global basis.
     */
    void evaluate(const Vec3&                      position_reference,
                  const Mat3&                      additional_rotation,
                  const VolumeStrainGreenLagrange& strain_global,
                  VolumeStressPK2&                 stress_global,
                  Mat6&                            tangent_global) const;

    /**
     * @brief Returns global tangent derivatives for additional rotation directions.
     */
    std::array<Mat6, 3> tangent_rotation_derivatives(
        const Vec3&                position_reference,
        const Mat3&                additional_rotation,
        const std::array<Mat3, 3>& additional_rotation_derivatives
    ) const;

    /**
     * @brief Outputs solid section details through the logger.
     */
    void info() override;

    /**
     * @brief Builds a compact one-line solid section summary.
     *
     * @return std::string Material, region and orientation.
     */
    std::string str() const override;

private:
    Mat3 section_orientation_basis(const Vec3& position_reference) const;
};
} // namespace fem
