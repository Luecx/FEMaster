/**
 * @file register_solid_section.inl
 * @brief Registers homogeneous Abaqus *SOLID SECTION definitions.
 *
 * The registration assigns homogeneous materials to continuum-solid or truss
 * element sets. Pure continuum sets map to `model::Model::solid_section()`, while
 * pure T3D2-derived truss sets map to `model::Model::truss_section()` and use the
 * section data value as cross-sectional area.
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
 * orientation is forwarded to continuum-solid sections. The section data value
 * is used as truss area for pure T3D2 sets. The target set must be non-empty and
 * contain only one supported section family.
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

        command.on_enter([material, elset, orientation](const fem::io::dsl::Keys& keys) {
            *material    = keys.raw("MATERIAL");
            *elset       = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("ATTRIBUTE").desc("Element-family section attribute")
                        .on_missing(fem::Precision{1}).on_empty(fem::Precision{1})
                )
                .bind([&model, material, elset, orientation](fem::Precision attribute) {
                    logging::error(model._data->elem_sets.has(*elset),
                        "SOLID SECTION element set '", *elset, "' does not exist");

                    auto region = model._data->elem_sets.get(*elset);
                    bool has_solid = false;
                    bool has_truss = false;
                    bool has_other = false;

                    for (fem::ID id : *region) {
                        auto& element = model._data->elements[id];
                        logging::error(element != nullptr,
                            "SOLID SECTION references undefined element ", id);

                        if (element->as<model::T3>() != nullptr) {
                            has_truss = true;
                            continue;
                        }

                        auto structural = element->as<model::StructuralElement>();
                        if (structural != nullptr && structural->is_solid()) {
                            has_solid = true;
                        } else {
                            has_other = true;
                        }
                    }

                    logging::error(region->size() > 0 && !has_other && !(has_solid && has_truss),
                        "SOLID SECTION requires a non-empty pure solid or pure truss element set");

                    if (has_truss) {
                        model.truss_section(*elset, *material, attribute);
                        return;
                    }

                    model.solid_section(*elset, *material, *orientation);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
