/**
 * @file register_solid_section.inl
 * @brief Registers solid section assignments.
 *
 * `SOLIDSECTION` binds a material to a target solid element region and may
 * attach an optional material coordinate system. The command resolves these
 * named model resources and constructs the `SolidSection` retained by semantic
 * Part topology.
 *
 * Constitutive evaluation, strain measures and stress recovery remain in the
 * solid elements and material implementation; this file only establishes their
 * section-level association.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../../model/model.h"
#include "../../../section/section_solid.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_solid_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SOLIDSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required()
                .key("ELSET").required()
                .key("ORIENTATION").optional()
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "SOLIDSECTION: no active part is available");

            const std::string material_name    = keys.raw("MATERIAL");
            const std::string elset_name       = keys.raw("ELSET");
            const std::string orientation_name = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};

            logging::error(part->elem_sets.has(elset_name),
                "SOLIDSECTION: element set ", elset_name, " is not defined in part ", part->name);
            logging::error(model._data->materials.has(material_name),
                "SOLIDSECTION: material ", material_name, " is not defined");
            logging::error(orientation_name.empty() || model._data->coordinate_systems.has(orientation_name),
                "SOLIDSECTION: coordinate system ", orientation_name, " is not defined");

            auto section = std::make_shared<SolidSection>();
            section->material_    = model._data->materials.get(material_name);
            section->region_      = part->elem_sets.get(elset_name);
            section->orientation_ = orientation_name.empty() ? nullptr : model._data->coordinate_systems.get(orientation_name);
            model.add_section(std::move(section));
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
