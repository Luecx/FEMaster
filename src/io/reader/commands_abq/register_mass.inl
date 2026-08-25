/**
 * @file register_mass.inl
 * @brief Registers Abaqus `*MASS` properties as `PointMassSection` assignments.
 *
 * Abaqus represents concentrated mass with two separate model concepts:
 * `*ELEMENT, TYPE=MASS` defines one-node element topology and `*MASS, ELSET=...`
 * assigns the physical mass property to a set of those elements. FEMaster maps
 * this directly to `model::PointElement` plus `PointMassSection`.
 *
 * The Abaqus reader replays the same input deck across multiple parser passes.
 * Part-local topology and section assignments are available before
 * `Model::compile()`, while assembly ELSETs are materialized only afterwards
 * because their references must be resolved through compiled Instance mappings.
 * Consequently the same `*MASS` record is encountered both before and after the
 * compilation boundary, but it must mutate the model in exactly one pass:
 *
 * - ROOT/PART `*MASS` is executed before compilation and stored as a part-local
 *   section assignment. `Model::compile()` later copies and remaps that section
 *   for every Instance of the Part.
 * - ASSEMBLY `*MASS` is only consumed before compilation, then executed in the
 *   post-compile assembly pass once the referenced assembly ELSET exists.
 * - On the opposite replay pass the already handled scope is still parsed, but
 *   its variant deliberately has no bind callback and therefore has no effect.
 *
 * Scope dispatch is expressed with `Condition::parent_is(...)` on the variants.
 * The model compile state is consulted only while constructing the registry to
 * decide which scope receives the active bind in the current replay. The bind
 * callbacks themselves contain only the operation for their semantic scope and
 * do not inspect parser stage or parent state.
 *
 * Only isotropic Abaqus mass is translated here. Anisotropic mass, orientation
 * and Abaqus damping parameters are outside the current implementation.
 *
 * @see model::PointElement
 * @see PointMassSection
 * @see model::Model::compile
 * @see commands::register_elset
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
 * Registers the Abaqus `*MASS, ELSET=...` command.
 *
 * Every accepted record contains one isotropic mass value and targets an ELSET
 * containing only `model::PointElement` objects. The command has separate DSL
 * variants for ROOT/PART and ASSEMBLY scope because those regions live in
 * different identifier spaces on opposite sides of `Model::compile()`.
 *
 * The registry is reconstructed for each reader replay. Before compilation the
 * ROOT/PART variant owns the bind and the ASSEMBLY variant is consume-only.
 * After compilation the roles are reversed: ROOT/PART is consume-only and the
 * ASSEMBLY variant owns the bind. This prevents duplicate section creation while
 * still allowing the complete deck to be parsed on every pass.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Model receiving part-local or compiled section assignments.
 */
inline void register_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("MASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Assign isotropic mass to Abaqus MASS elements in an ELSET.");

        // The keyword is replayed in multiple parser passes. Keep its ELSET name
        // in command-local state so whichever variant is active for this replay
        // can resolve the target in the appropriate identifier space.
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

        // ------------------------------------------------------------------
        // Pre-compile replay
        // ------------------------------------------------------------------
        //
        // Part-local ELSETs and PointElements already exist at this point, so
        // ROOT/PART MASS properties can be stored immediately as semantic
        // section assignments on the active Part. Assembly ELSETs do not yet
        // exist in compiled identifier space and therefore cannot be resolved.
        if (!model._data->compiled) {
            // ROOT/PART: execute the MASS property now.
            //
            // The region belongs to the active semantic Part. Model::add_section()
            // therefore stores the PointMassSection on that Part, and compile()
            // later copies/remaps it once for every Instance.
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

                        // `*MASS` is a property of Abaqus MASS elements. Reject a
                        // mixed ELSET instead of silently assigning the section to
                        // unrelated structural element formulations.
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

            // ASSEMBLY: consume only during the pre-compile replay.
            //
            // The parser must still recognize the data line so the deck remains
            // syntactically valid, but no bind is attached because assembly
            // ELSETs are materialized only after Model::compile(). The same
            // record will be replayed later and handled by the active ASSEMBLY
            // variant below.
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

        // ------------------------------------------------------------------
        // Post-compile assembly replay
        // ------------------------------------------------------------------
        //
        // Part-local MASS sections have already been expanded by compile(). The
        // compiled assembly ELSET namespace is now available, so this replay must
        // ignore ROOT/PART records and execute only ASSEMBLY records.

        // ROOT/PART: consume only after compilation.
        //
        // These records already created their semantic PointMassSection during
        // the pre-compile replay. Replaying their bind here would create a second
        // property assignment, so the variant intentionally has no callback.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
            )
        );

        // ASSEMBLY: execute the MASS property after compilation.
        //
        // The target ELSET now contains dense compiled element ids. The section
        // is therefore created directly in assembly space and Model::add_section()
        // assigns it immediately to the referenced compiled PointElements.
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

                    // Apply the same topology validation as in semantic Part
                    // space, now against dense compiled element storage.
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
