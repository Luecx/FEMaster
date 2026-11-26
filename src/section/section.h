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
#include "../core/types_eig.h"
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

    // Positive section id. It is created once per section definition and kept
    // by instance copies during model compilation, e.g. for the FRD group field.
    Index id_ = create_id();

    // Creates a new unique positive section id, starting at 1
    static Index create_id() {
        static Index next_id = 0;
        return ++next_id;
    }

    /**
     * @brief Defaulted virtual destructor for polymorphic deletion.
     */
    virtual ~Section() = default;

    // Material evaluated at the element-local material point mp of every
    // integration point. Single-material sections use material_ everywhere.
    virtual material::Material::Ptr material_at(Index mp) const {
        (void) mp;
        return material_;
    }

    // Largest constitutive state width required by any material point of this
    // section. Sizes the model-wide ELEMENT_MP state field.
    virtual Index material_state_size() const {
        return material_ && material_->has_elasticity() ? material_->elasticity()->state_size() : 0;
    }

    // Independent copy of this section for one compiled instance. The copy keeps
    // the section id, targets the given compiled element region and places all
    // spatial section data rigidly according to x = R X + t.
    virtual Ptr instantiate(model::ElementRegion::Ptr region,
                            const Mat3&               rotation,
                            const Vec3&               translation) const = 0;

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
