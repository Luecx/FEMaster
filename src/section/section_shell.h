/**
 * @file section_shell.h
 * @brief Declares shell section properties and material matrices.
 *
 * @see src/section/section_shell.cpp
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "section.h"

#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"

namespace fem {
/**
 * @struct ShellSection
 * @brief Section definition for shell elements.
 */
struct ShellSection : Section {
    using Ptr = std::shared_ptr<ShellSection>; ///< Shared pointer alias for shell sections.

    Precision                  thickness_   = 1.0;     ///< Shell thickness.
    cos::CoordinateSystem::Ptr orientation_ = nullptr; ///< Optional section/material orientation.

    /**
     * @brief Computes the shell ABD matrix from material elasticity and thickness.
     *
     * @return StaticMatrix<6, 6> Membrane-bending stiffness matrix.
     */
    virtual StaticMatrix<6, 6> get_abd();

    /**
     * @brief Computes the shell transverse shear matrix.
     *
     * @return StaticMatrix<2, 2> Shear stiffness matrix.
     */
    virtual StaticMatrix<2, 2> get_shear();

    /**
     * @brief Checks whether the section material defines density.
     *
     * @return bool True when material density is available.
     */
    bool has_density() const;

    /**
     * @brief Returns material density or zero when missing.
     *
     * @return Precision Density value.
     */
    Precision get_density() const;

    /**
     * @brief Outputs shell section details through the logger.
     */
    void info() override;

    /**
     * @brief Builds a compact one-line shell section summary.
     *
     * @return std::string Material, region, density, orientation and thickness.
     */
    std::string str() const override;
};
} // namespace fem
