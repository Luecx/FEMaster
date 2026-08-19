/**
 * @file register_shell_section.inl
 * @brief Registers homogeneous Abaqus shell sections.
 *
 * `SHELL SECTION` resolves the target element set, material, thickness and
 * optional material orientation and translates them into an integrated
 * FEMaster shell section. Section data remains associated with the active Part
 * until model compilation expands its Instances.
 *
 * Layered composites and unsupported Abaqus section controls are rejected or
 * deliberately left outside this homogeneous-section mapping.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>

#include "../../../model/model.h"
#include "../../../section/section_shell_integrated.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_shell_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SHELLSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));

        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").required()
                .key("MATERIAL").required()
                .key("ORIENTATION").optional()
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
                    .one<Precision>().name("THICKNESS")
                    .one<int>().name("INTEGRATIONPOINTS").on_missing(5).on_empty(5)
                )
                .bind([&model, material, elset, orientation](Precision thickness, int integration_points) {
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SHELL SECTION: no active part is available");
                    logging::error(integration_points == 5,
                        "SHELL SECTION: exactly 5 Simpson integration points are supported");
                    logging::error(thickness > Precision(0),
                        "SHELL SECTION: thickness must be positive");
                    logging::error(part->elem_sets.has(*elset),
                        "SHELL SECTION: element set ", *elset, " does not exist in part ", part->name);
                    logging::error(model._data->materials.has(*material),
                        "SHELL SECTION: material ", *material, " does not exist");
                    logging::error(orientation->empty() || model._data->coordinate_systems.has(*orientation),
                        "SHELL SECTION: orientation ", *orientation, " does not exist");

                    model.add_section(std::make_shared<IntegratedShellSection>(
                        model._data->materials.get(*material),
                        part->elem_sets.get(*elset),
                        thickness,
                        orientation->empty() ? nullptr : model._data->coordinate_systems.get(*orientation),
                        0
                    ));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
