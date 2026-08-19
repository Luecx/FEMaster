/**
 * @file register_solid_section.inl
 * @brief Registers Abaqus solid and truss section definitions.
 *
 * `SOLID SECTION` binds a material and optional orientation to a solid element
 * set, while the supported truss section form additionally reads its area. The
 * command creates the corresponding FEMaster section object on the active
 * semantic Part.
 *
 * Constitutive integration and element kinematics remain in the concrete
 * sections and elements; this file owns only Abaqus syntax translation and
 * resource resolution.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>
#include <string>

#include "../../../model/model.h"
#include "../../../model/truss/truss.h"
#include "../../../section/section_solid.h"
#include "../../../section/section_truss.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_solid_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SOLIDSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));

        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto has_truss   = std::make_shared<bool>(false);
        auto truss_created = std::make_shared<bool>(false);

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").required()
                .key("MATERIAL").required()
                .key("ORIENTATION").optional()
        );
        command.on_enter([&model, material, elset, orientation, has_truss, truss_created](const fem::io::dsl::Keys& keys) {
            *material      = keys.raw("MATERIAL");
            *elset         = keys.raw("ELSET");
            *orientation   = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *has_truss     = false;
            *truss_created = false;

            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "SOLID SECTION: no active part is available");
            logging::error(part->elem_sets.has(*elset),
                "SOLID SECTION: element set ", *elset, " does not exist in part ", part->name);
            logging::error(model._data->materials.has(*material),
                "SOLID SECTION: material ", *material, " does not exist");
            logging::error(orientation->empty() || model._data->coordinate_systems.has(*orientation),
                "SOLID SECTION: orientation ", *orientation, " does not exist");

            const auto region = part->elem_sets.get(*elset);
            bool has_solid = false;
            bool has_other = false;
            for (const ID id : *region) {
                const auto it = part->elements.find(id);
                logging::error(it != part->elements.end() && it->second != nullptr,
                    "SOLID SECTION: element ", id, " does not exist in part ", part->name);

                if (it->second->as<model::T3>()) {
                    *has_truss = true;
                } else if (auto* structural = it->second->as<model::StructuralElement>(); structural && structural->is_solid()) {
                    has_solid = true;
                } else {
                    has_other = true;
                }
            }

            logging::error(region->size() > 0 && !has_other && !(has_solid && *has_truss),
                "SOLID SECTION: element set must be a non-empty pure solid or pure truss set");

            if (has_solid) {
                auto section = std::make_shared<SolidSection>();
                section->material_    = model._data->materials.get(*material);
                section->region_      = region;
                section->orientation_ = orientation->empty() ? nullptr : model._data->coordinate_systems.get(*orientation);
                model.add_section(std::move(section));
            }
        });
        command.on_exit([has_truss, truss_created](const fem::io::dsl::Keys&) {
            logging::error(!*has_truss || *truss_created,
                "SOLID SECTION: truss section requires a cross-sectional area");
        });
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("ATTRIBUTE").on_missing(Precision{1}).on_empty(Precision{1})
                )
                .bind([&model, material, elset, has_truss, truss_created](Precision area) {
                    if (!*has_truss) return;
                    logging::error(area > Precision(0),
                        "SOLID SECTION: truss area must be positive");
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SOLID SECTION: no active part is available");
                    model.add_section(std::make_shared<TrussSection>(
                        model._data->materials.get(*material),
                        part->elem_sets.get(*elset),
                        area
                    ));
                    *truss_created = true;
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
