/**
 * @file register_spring.inl
 * @brief Registers linear SPRING1 point-element properties.
 *
 * `SPRING, ELSET=...` assigns one scalar ground stiffness to one-node `SPRING1`
 * topology. Degrees of freedom 1--3 map to translational stiffness and 4--6 to
 * rotational stiffness in `PointMassSection`.
 *
 * ROOT/PART records resolve semantic Part-local ELSETs before
 * `Model::compile()`, while ASSEMBLY records resolve dense compiled ELSETs
 * afterwards. Both scope variants are registered once; the parsed-deck reader
 * controls when each occurrence executes.
 *
 * @see model::PointElement
 * @see PointMassSection
 * @see commands::register_point_mass
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "../../../model/element/point.h"
#include "../../../model/model.h"
#include "../../../section/section_point_mass.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <memory>
#include <string>

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers constant linear SPRING1 properties in Part and assembly space.
 *
 * The two variants share the Abaqus two-line data layout. Explicit compile-state
 * checks document the required semantic processing order and reject accidental
 * execution on the wrong side of the `Model::compile()` boundary.
 */
inline void register_spring(dsl::Registry& registry, model::Model& model) {
    registry.command("SPRING", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Assign constant linear ground stiffness to SPRING1 point elements; native legacy alternative: POINTMASS.");

        auto elset = std::make_shared<std::string>();

        command.keyword(
            dsl::KeywordSpec::make()
                .key("ELSET").required().doc("Element set containing SPRING1 point elements")
        );

        command.on_enter([elset](const dsl::Keys& keys) {
            *elset = keys.raw("ELSET");
        });

        // Part-local stiffness is attached before compile() and copied through
        // every Instance by the ordinary section-expansion path.
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(2).max(2))
                .pattern(dsl::Pattern::make()
                    .allow_multiline()
                    .one<ID>().name("DOF").desc("Nodal degree of freedom 1..6")
                    .one<Precision>().name("STIFFNESS").desc("Linear spring stiffness")
                )
                .bind([&model, elset](ID dof, Precision stiffness) {
                    logging::error(!model._data->compiled,
                        "SPRING: ROOT/PART assignment must be processed before model compilation");

                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SPRING: no active part is available");
                    logging::error(part->elem_sets.has(*elset),
                        "SPRING: element set ", *elset,
                        " does not exist in part ", part->name);

                    const auto region = part->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "SPRING: element set ", *elset,
                        " is empty in part ", part->name);

                    for (const ID id : *region) {
                        const auto it = part->elements.find(id);
                        logging::error(it != part->elements.end() && it->second != nullptr,
                            "SPRING: element ", id,
                            " does not exist in part ", part->name);
                        logging::error(it->second->as<model::PointElement>() != nullptr,
                            "SPRING: element set ", *elset,
                            " contains non-point element ", id);
                    }

                    logging::error(dof >= 1 && dof <= 6,
                        "SPRING: degree of freedom must be between 1 and 6");

                    Vec3 spring        = Vec3::Zero();
                    Vec3 rotary_spring = Vec3::Zero();
                    if (dof <= 3) {
                        spring(static_cast<Index>(dof - 1)) = stiffness;
                    } else {
                        rotary_spring(static_cast<Index>(dof - 4)) = stiffness;
                    }

                    model.add_section(std::make_shared<PointMassSection>(
                        region, Precision(0), Vec3::Zero(), spring, rotary_spring));
                })
            )
        );

        // Assembly stiffness resolves dense compiled element ids after compile().
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::parent_is("ASSEMBLY"))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(2).max(2))
                .pattern(dsl::Pattern::make()
                    .allow_multiline()
                    .one<ID>().name("DOF").desc("Nodal degree of freedom 1..6")
                    .one<Precision>().name("STIFFNESS").desc("Linear spring stiffness")
                )
                .bind([&model, elset](ID dof, Precision stiffness) {
                    logging::error(model._data->compiled,
                        "SPRING: ASSEMBLY assignment requires a compiled model");
                    logging::error(model._data->elem_sets.has(*elset),
                        "SPRING: compiled element set ", *elset, " does not exist");

                    const auto region = model._data->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "SPRING: compiled element set ", *elset, " is empty");

                    for (const ID id : *region) {
                        logging::error(id >= 0
                                    && static_cast<std::size_t>(id) < model._data->elements.size()
                                    && model._data->elements[static_cast<std::size_t>(id)] != nullptr,
                            "SPRING: compiled element ", id, " does not exist");
                        logging::error(model._data->elements[static_cast<std::size_t>(id)]->as<model::PointElement>() != nullptr,
                            "SPRING: compiled element set ", *elset,
                            " contains non-point element ", id);
                    }

                    logging::error(dof >= 1 && dof <= 6,
                        "SPRING: degree of freedom must be between 1 and 6");

                    Vec3 spring        = Vec3::Zero();
                    Vec3 rotary_spring = Vec3::Zero();
                    if (dof <= 3) {
                        spring(static_cast<Index>(dof - 1)) = stiffness;
                    } else {
                        rotary_spring(static_cast<Index>(dof - 4)) = stiffness;
                    }

                    model.add_section(std::make_shared<PointMassSection>(
                        region, Precision(0), Vec3::Zero(), spring, rotary_spring));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
