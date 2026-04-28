/**
 * @file section_truss.h
 * @brief Declares truss section properties.
 *
 * @see src/section/section_truss.cpp
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "section.h"

namespace fem {
/**
 * @struct TrussSection
 * @brief Section definition for truss elements.
 */
struct TrussSection : Section {
    using Ptr = std::shared_ptr<TrussSection>; ///< Shared pointer alias for truss sections.

    Precision area_ = Precision(0); ///< Cross-section area.

    /**
     * @brief Outputs truss section details through the logger.
     */
    void info() override;

    /**
     * @brief Builds a compact one-line truss section summary.
     *
     * @return std::string Material, region and area.
     */
    std::string str() const override;
};
} // namespace fem
