/**
 * @file section_beam.h
 * @brief Declares beam section properties.
 *
 * @see src/section/section_beam.cpp
 * @see src/section/profile.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "section.h"
#include "profile.h"

#include "../core/types_eig.h"

namespace fem {
/**
 * @struct BeamSection
 * @brief Section definition for beam elements.
 */
struct BeamSection : Section {
    using Ptr = std::shared_ptr<BeamSection>; ///< Shared pointer alias for beam sections.

    Vec3         direction_ = Vec3::Zero(); ///< Optional beam section direction vector.
    Profile::Ptr profile_  = nullptr;       ///< Profile associated with the beam section.

    /**
     * @brief Outputs beam section details through the logger.
     */
    void info() override;

    /**
     * @brief Builds a compact one-line beam section summary.
     *
     * @return std::string Material, region, profile and direction.
     */
    std::string str() const override;
};
} // namespace fem
