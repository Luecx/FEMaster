/**
 * @file register_solid_section.inl
 * @brief Registers homogeneous Abaqus *SOLID SECTION definitions.
 *
 * The registration assigns homogeneous materials to continuum-solid or truss
 * element sets. Pure continuum sets map directly to
 * `model::Model::solid_section()` without requiring a data line. Pure
 * T3D2-derived truss sets consume one data value as their cross-sectional area
 * and map to `model::Model::truss_section()`.
 *
 * @see model::Model::solid_section
 * @see model::Model::truss_section
 * @see model::T3
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <memory>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../../model/truss/truss.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers homogeneous Abaqus solid and truss section assignment.
 *
 * `ELSET` and `MATERIAL` identify the target elements and material. An optional
 * orientation is forwarded to continuum-solid sections, which are created when
 * the keyword is entered. Pure T3D2 sets defer section creation until their
 * cross-sectional area is read from the conditionally required command data
 * line. The target set must be non-empty and contain only one supported section
 * family.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the section definition.
 */
inline void register_solid_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SOLIDSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a homogeneous Abaqus solid or truss section to an element set.");

        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();

        auto has_truss             = std::make_shared<bool>(false);
        auto truss_section_created = std::make_shared<bool>(false);

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET")
                    .required()
                    .doc("Target element set")
                .key("MATERIAL")
                    .required()
                    .doc("Material assigned to the section")
                .key("ORIENTATION")
                    .optional()
                    .doc("Optional continuum-solid material orientation")
        );

        command.on_enter([&model, material, elset, orientation, has_truss, truss_section_created](
            const fem::io::dsl::Keys& keys
        ) {
            // Store the section definition and reset the data-line state
            *material    = keys.raw("MATERIAL");
            *elset       = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};

            *has_truss             = false;
            *truss_section_created = false;

            // Classify the target set before creating an element-family section
            logging::error(model._data->elem_sets.has(*elset),
                "SOLID SECTION element set '", *elset, "' does not exist");

            auto region = model._data->elem_sets.get(*elset);
            bool has_solid = false;
            bool has_other = false;

            for (fem::ID id : *region) {
                auto& element = model._data->elements[id];
                logging::error(element != nullptr,
                    "SOLID SECTION references undefined element ", id);

                if (element->as<model::T3>() != nullptr) {
                    *has_truss = true;
                    continue;
                }

                auto structural = element->as<model::StructuralElement>();
                if (structural != nullptr && structural->is_solid()) {
                    has_solid = true;
                } else {
                    has_other = true;
                }
            }

            logging::error(region->size() > 0 && !has_other && !(has_solid && *has_truss),
                "SOLID SECTION requires a non-empty pure solid or pure truss element set");

            // Continuum solids require no section data and can be assigned immediately
            if (has_solid) {
                model.solid_section(*elset, *material, *orientation);
            }
        });

        command.on_exit([has_truss, truss_section_created](const fem::io::dsl::Keys&) {
            // A truss section is incomplete until its cross-sectional area was read
            logging::error(!*has_truss || *truss_section_created,
                "SOLID SECTION for a truss element set requires a cross-sectional area");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("ATTRIBUTE").desc("Element-family section attribute")
                        .on_missing(fem::Precision{1}).on_empty(fem::Precision{1})
                )
                .bind([&model, material, elset, has_truss, truss_section_created](
                    fem::Precision attribute
                ) {
                    // Interpret section data only as the area of a pure truss set
                    if (*has_truss) {
                        model.truss_section(*elset, *material, attribute);
                        *truss_section_created = true;
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
