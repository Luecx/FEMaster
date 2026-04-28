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

namespace fem {
/**
 * @struct SolidSection
 * @brief Section definition for solid elements.
 */
struct SolidSection : Section {
    using Ptr = std::shared_ptr<SolidSection>; ///< Shared pointer alias for solid sections.

    cos::CoordinateSystem::Ptr orientation_ = nullptr; ///< Optional material orientation.

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
};
} // namespace fem
