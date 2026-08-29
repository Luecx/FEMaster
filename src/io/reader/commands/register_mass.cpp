/**
 * @file register_mass.cpp
 * @brief Registers shared `*MASS` point-element properties.
 *
 * `*ELEMENT, TYPE=MASS` defines one-node point topology while
 * `*MASS, ELSET=...` assigns its isotropic translational mass through
 * `PointMassSection`. ROOT/PART records operate in semantic Part-local id space;
 * ASSEMBLY records operate in dense compiled assembly space.
 *
 * Both scope variants are registered once. The parsed-deck reader controls when
 * each occurrence executes: Part-local assignments before `Model::compile()` and
 * assembly assignments afterwards. The command therefore contains no parser
 * stage or registration-time compile-state branching.
 *
 * @see model::PointElement
 * @see PointMassSection
 * @see model::Model::compile
 * @see commands::register_elset
 * @see commands::register_point_mass
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../model/element/point.h"
#include "../../../model/model.h"
#include "../../../section/section_point_mass.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <memory>
#include <string>

namespace fem::io::reader::commands {

/**
 * @brief Registers isotropic mass assignments for Part-local and assembly ELSETs.
 *
 * The two variants share one syntax but resolve their target in different model
 * namespaces. Explicit compile-state checks document the required semantic
 * processing order and reject accidental execution on the wrong side of the
 * `Model::compile()` boundary.
 */
void register_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("MASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Assign isotropic mass to MASS point elements in an ELSET; native legacy alternative: POINTMASS.");

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

        // Part-local properties are attached before compile() and copied through
        // every Instance by the ordinary section-expansion path.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
                .bind([&model, elset](Precision mass) {
                    logging::error(!model._data->compiled,
                        "MASS: ROOT/PART assignment must be processed before model compilation");

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

        // Assembly properties resolve already materialized dense ELSET ids after
        // compile() and are assigned directly to the compiled point elements.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is("ASSEMBLY"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
                .bind([&model, elset](Precision mass) {
                    logging::error(model._data->compiled,
                        "MASS: ASSEMBLY assignment requires a compiled model");
                    logging::error(model._data->elem_sets.has(*elset),
                        "MASS: compiled element set ", *elset, " does not exist");

                    const auto region = model._data->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "MASS: compiled element set ", *elset, " is empty");

                    for (const ID id : *region) {
                        logging::error(id >= 0
                                    && static_cast<std::size_t>(id) < model._data->elements.size()
                                    && model._data->elements[static_cast<std::size_t>(id)] != nullptr,
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

} // namespace fem::io::reader::commands
