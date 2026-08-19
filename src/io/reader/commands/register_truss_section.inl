/**
 * @file register_truss_section.inl
 * @brief Registers truss section assignments.
 *
 * `TRUSSSECTION` resolves a material and target truss element region, reads the
 * positive cross-sectional area and creates the corresponding `TrussSection`.
 * The association is retained on semantic topology before dense assembly data
 * is generated.
 *
 * Axial kinematics, constitutive response and area scaling are evaluated later
 * by the truss element and section during analysis.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>

#include "../../../model/model.h"
#include "../../../section/section_truss.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_truss_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("TRUSSSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));

        auto material = std::make_shared<std::string>();
        auto elset    = std::make_shared<std::string>();
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required()
                .key("ELSET").required()
        );
        command.on_enter([material, elset](const fem::io::dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset    = keys.raw("ELSET");
        });
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("AREA")
                )
                .bind([&model, material, elset](Precision area) {
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "TRUSSSECTION: no active part is available");
                    logging::error(part->elem_sets.has(*elset),
                        "TRUSSSECTION: element set ", *elset, " is not defined in part ", part->name);
                    logging::error(model._data->materials.has(*material),
                        "TRUSSSECTION: material ", *material, " is not defined");
                    logging::error(area > Precision(0),
                        "TRUSSSECTION: area must be positive");

                    model.add_section(std::make_shared<TrussSection>(
                        model._data->materials.get(*material),
                        part->elem_sets.get(*elset),
                        area
                    ));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
