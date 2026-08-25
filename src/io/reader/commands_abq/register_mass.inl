/**
 * @file register_mass.inl
 * @brief Registers isotropic Abaqus MASS properties as point-mass sections.
 *
 * Abaqus separates one-node `MASS` topology from its scalar property. FEMaster
 * mirrors that model directly: `model::PointElement` owns element identity and
 * connectivity, while `PointMassSection` supplies the concentrated mass.
 *
 * Scope selection is expressed in the DSL variants themselves. During the
 * pre-compile parser pass only ROOT/PART variants mutate semantic sections;
 * during the post-compile pass only the ASSEMBLY variant mutates compiled
 * sections. The bind callbacks therefore contain no parser-stage branching.
 *
 * Only the standard isotropic Abaqus form is supported here. Anisotropic mass,
 * orientation and damping parameters remain outside this translation.
 *
 * @see model::PointElement
 * @see PointMassSection
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include <memory>
#include <string>

#include "../../../model/element/point.h"
#include "../../../model/model.h"
#include "../../../section/section_point_mass.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers `*MASS, ELSET=...` as a section assignment on PointElements.
 *
 * The registry is rebuilt for every parser pass. The model compile state is
 * therefore used only while constructing the DSL: the active semantic pass gets
 * the ROOT/PART bind and a consuming ASSEMBLY variant, while the compiled pass
 * gets the inverse. Runtime bind callbacks select behavior solely through their
 * `parent_is(...)` variant.
 *
 * @param registry Parser registry receiving the MASS command.
 * @param model Model receiving semantic or compiled section assignments.
 */
inline void register_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("MASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Assign isotropic mass to Abaqus MASS elements in an ELSET.");

        auto elset = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").required().doc("Element set containing MASS elements")
                .key("TYPE").optional("ISOTROPIC").allowed({"ISOTROPIC"})
                    .doc("Only isotropic point mass is currently supported")
        );

        command.on_enter([elset](const fem::io::dsl::Keys& keys) {
            *elset = keys.raw("ELSET");
        });

        if (!model._data->compiled) {
            command.variant(fem::io::dsl::Variant::make()
                .when(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}))
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(1).max(1))
                    .pattern(fem::io::dsl::Pattern::make()
                        .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                    )
                    .bind([&model, elset](Precision mass) {
                        const auto part = model._data->parts.get();
                        logging::error(part != nullptr,
                            "MASS: no active part is available");
                        logging::error(part->elem_sets.has(*elset),
                            "MASS: element set ", *elset, " does not exist in part ", part->name);

                        const auto region = part->elem_sets.get(*elset);
                        logging::error(region != nullptr && region->size() > 0,
                            "MASS: element set ", *elset, " is empty in part ", part->name);

                        for (const ID id : *region) {
                            const auto it = part->elements.find(id);
                            logging::error(it != part->elements.end() && it->second != nullptr,
                                "MASS: element ", id, " does not exist in part ", part->name);
                            logging::error(it->second->as<model::PointElement>() != nullptr,
                                "MASS: element set ", *elset, " contains non-MASS element ", id);
                        }

                        model.add_section(std::make_shared<PointMassSection>(region, mass));
                    })
                )
            );

            command.variant(fem::io::dsl::Variant::make()
                .when(fem::io::dsl::Condition::parent_is("ASSEMBLY"))
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(1).max(1))
                    .pattern(fem::io::dsl::Pattern::make()
                        .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                    )
                )
            );
            return;
        }

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is("ASSEMBLY"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
                .bind([&model, elset](Precision mass) {
                    logging::error(model._data->elem_sets.has(*elset),
                        "MASS: compiled element set ", *elset, " does not exist");

                    const auto region = model._data->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "MASS: compiled element set ", *elset, " is empty");

                    for (const ID id : *region) {
                        logging::error(id >= 0 && static_cast<std::size_t>(id) < model._data->elements.size() && model._data->elements[static_cast<std::size_t>(id)] != nullptr,
                            "MASS: compiled element ", id, " does not exist");
                        logging::error(model._data->elements[static_cast<std::size_t>(id)]->as<model::PointElement>() != nullptr,
                            "MASS: compiled element set ", *elset, " contains non-MASS element ", id);
                    }

                    model.add_section(std::make_shared<PointMassSection>(region, mass));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
