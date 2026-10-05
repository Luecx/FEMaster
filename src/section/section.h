/**
 * @file section.h
 * @brief Declares the base interface for FEM section definitions.
 *
 * Sections bind element regions to material data and provide type-specific
 * section properties for solids, shells, beams and trusses.
 *
 * @see src/section/section.cpp
 * @see src/section/section_solid.h
 * @see src/section/section_shell.h
 * @see src/section/section_beam.h
 * @see src/section/section_truss.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "../core/castable.h"
#include "../core/printable.h"
#include "../data/region.h"
#include "../material/material.h"

#include <memory>
#include <string>

namespace fem {
/**
 * @struct Section
 * @brief Abstract base class for FEM sections.
 *
 * Every concrete section owns the material-region pairing used by elements and
 * adds only the extra properties needed by that element family. Castable
 * provides mutable and const runtime access to concrete section types.
 */
struct Section : public fem::Printable, public Castable {
    using Ptr = std::shared_ptr<Section>; ///< Shared pointer alias for section ownership.

    material::Material::Ptr  material_ = nullptr; ///< Material associated with the section.
    model::ElementRegion::Ptr region_  = nullptr; ///< Element region associated with the section.

    /**
     * @brief Defaulted virtual destructor for polymorphic deletion.
     */
    virtual ~Section() = default;

    /**
     * @brief Outputs section details through the logger.
     */
    virtual void info() = 0;

    /**
     * @brief Builds a compact one-line section summary.
     *
     * @return std::string Section type, material, region and properties.
     */
    std::string str() const override = 0;
};
} // namespace fem
