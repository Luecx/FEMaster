/**
 * @file register_mass.inl
 * @brief Registers the Abaqus `*MASS` property command in the replay-based reader.
 *
 * This file translates Abaqus concentrated-mass properties into FEMaster section
 * assignments. Abaqus separates the one-node `MASS` element topology from the
 * property assigned to that topology: `*ELEMENT, TYPE=MASS` defines element ids,
 * connectivity and ELSET membership, while `*MASS, ELSET=...` assigns the scalar
 * isotropic mass. FEMaster represents the same separation with
 * `model::PointElement` and `PointMassSection`.
 *
 * The Abaqus reader processes the same input deck in multiple semantic passes
 * around the one-way `Model::compile()` boundary. This is relevant for `*MASS`
 * because the target ELSET lives in different identifier spaces depending on
 * its parent scope:
 *
 * - ROOT/PART records reference semantic Part-local element sets. They can be
 *   resolved before compilation and are stored as Part-local section
 *   assignments. `Model::compile()` later copies and remaps those assignments
 *   for every Instance of the Part.
 * - ASSEMBLY records reference assembly element sets. Those sets are
 *   materialized only after compilation, when local Instance element ids can be
 *   resolved to dense assembly ids.
 *
 * The command therefore uses two DSL variants in each relevant reader replay.
 * One variant performs the model mutation for the scope that is resolvable in
 * that replay, while the other variant only consumes and validates the same
 * input syntax. This prevents duplicate section creation although the complete
 * deck is parsed more than once.
 *
 * Parent-scope dispatch belongs to the DSL and is expressed exclusively through
 * `Condition::parent_is(...)` on each variant. The bind callbacks themselves do
 * not inspect parser parent state. The model compile state is used only while
 * registering the command to select which replay receives the mutating bind.
 *
 * This translation currently supports only the isotropic Abaqus `*MASS` form.
 * Anisotropic mass, orientation-dependent mass and Abaqus damping parameters are
 * not represented here. Mass, rotary inertia and ground-spring behavior of the
 * FEMaster `PointMassSection` remain responsibilities of the section and point
 * element implementations rather than of this parser command.
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
 * Registers `*MASS, ELSET=...` for Part-local and assembly-level MASS elements.
 *
 * The command accepts one isotropic mass value and requires an ELSET containing
 * only `model::PointElement` objects. The target region is validated before a
 * `PointMassSection` is created so a mixed ELSET cannot silently assign point
 * properties to unrelated structural formulations.
 *
 * Registration depends on which side of `Model::compile()` the current parser
 * replay runs:
 *
 * 1. Before compilation, ROOT/PART records are executed against the active
 *    semantic Part. ASSEMBLY records are parsed without a bind because their
 *    assembly ELSETs do not exist yet.
 * 2. After compilation, ROOT/PART records are parsed without a bind because
 *    their sections have already been copied into assembly space. ASSEMBLY
 *    records are now executed against compiled ELSETs and dense element ids.
 *
 * The consume-only variants are necessary because the Abaqus reader replays the
 * complete input deck. They preserve syntax handling in every pass while making
 * the section assignment itself occur exactly once. Scope routing remains
 * declarative through `parent_is(...)`; no bind callback branches on parser
 * stage or parent scope.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Model receiving either semantic Part-local or compiled assembly
 *              section assignments.
 */
inline void register_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("MASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Assign isotropic mass to Abaqus MASS elements in an ELSET.");

        // The ELSET name belongs to one keyword occurrence but is needed later
        // when the selected data variant executes. Shared command-local storage
        // therefore carries it from keyword parsing into the variant bind.
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
        // ROOT/PART topology and ELSETs already exist in their semantic Part
        // identifier spaces. Assembly ELSETs, in contrast, are intentionally
        // deferred until after compile() because Instance-local element ids must
        // first be mapped into the dense assembly namespace.
        if (!model._data->compiled) {
            // ROOT/PART variant -- execute before compile().
            //
            // This is the mutating variant for semantic input. `parent_is(...)`
            // guarantees that the bind is entered only for a MASS record whose
            // parent is ROOT or PART; the callback therefore does not need to
            // inspect any parser stage or scope flag. The referenced ELSET is
            // resolved directly from the active Part and validated to contain
            // only PointElements. `Model::add_section()` sees an uncompiled model
            // and stores the resulting PointMassSection on that Part. During
            // compile(), the ordinary section-copy path creates one remapped
            // assignment for every Instance of the Part.
            command.variant(fem::io::dsl::Variant::make()
                .when(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}))
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(1).max(1))
                    .pattern(fem::io::dsl::Pattern::make()
                        .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                    )
                    .bind([&model, elset](Precision mass) {
                        // Resolve the semantic target region from the currently
                        // active Part. At this stage ids are still Part-local.
                        const auto part = model._data->parts.get();
                        logging::error(part != nullptr,
                            "MASS: no active part is available");
                        logging::error(part->elem_sets.has(*elset),
                            "MASS: element set ", *elset, " does not exist in part ", part->name);

                        const auto region = part->elem_sets.get(*elset);
                        logging::error(region != nullptr && region->size() > 0,
                            "MASS: element set ", *elset, " is empty in part ", part->name);

                        // Abaqus *MASS is a property of MASS elements. Validate
                        // the complete region before creating the section so a
                        // partially valid mixed ELSET never mutates the model.
                        for (const ID id : *region) {
                            const auto it = part->elements.find(id);
                            logging::error(it != part->elements.end() && it->second != nullptr,
                                "MASS: element ", id, " does not exist in part ", part->name);
                            logging::error(it->second->as<model::PointElement>() != nullptr,
                                "MASS: element set ", *elset, " contains non-MASS element ", id);
                        }

                        // add_section() stores this assignment on the active Part
                        // because the model has not crossed compile() yet.
                        model.add_section(std::make_shared<PointMassSection>(region, mass));
                    })
                )
            );

            // ASSEMBLY variant -- consume only before compile().
            //
            // The complete deck is already being replayed, so an assembly-level
            // MASS record must still match its normal data grammar. Its target
            // ELSET cannot be resolved yet, however, because assembly ELSETs are
            // built only after Instance element mappings exist. This variant
            // therefore deliberately has no bind callback: it consumes the data
            // line without modifying the model. The same record is encountered
            // again in the post-compile assembly replay, where the corresponding
            // ASSEMBLY variant below performs the actual section assignment.
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
        // Part-local PointMassSections have already been expanded and assigned
        // to their compiled Instance elements by compile(). Assembly ELSETs now
        // exist in dense identifier space, so the mutating responsibility moves
        // from ROOT/PART to ASSEMBLY for this replay.

        // ROOT/PART variant -- consume only after compile().
        //
        // These records were already executed in the pre-compile replay and the
        // resulting Part-local sections were copied by compile(). Because the
        // complete input deck is replayed again, the parser still needs a valid
        // ROOT/PART grammar here, but attaching the original bind would create a
        // duplicate section assignment. The missing bind is therefore deliberate:
        // this variant recognizes and consumes the record without touching the
        // compiled model.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
            )
        );

        // ASSEMBLY variant -- execute after compile().
        //
        // Assembly ELSETs have now been materialized from their Instance/local
        // references into dense global element ids. `parent_is("ASSEMBLY")`
        // routes only assembly-level MASS records into this bind. The callback
        // validates that every referenced compiled element exists and is a
        // PointElement, then creates the PointMassSection directly in assembly
        // space. Because the model is compiled, `Model::add_section()` stores and
        // immediately assigns the section to the referenced compiled elements.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is("ASSEMBLY"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
                .bind([&model, elset](Precision mass) {
                    // Resolve and validate the target against the compiled
                    // assembly namespace before creating any section object.
                    logging::error(model._data->elem_sets.has(*elset),
                        "MASS: compiled element set ", *elset, " does not exist");

                    const auto region = model._data->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "MASS: compiled element set ", *elset, " is empty");

                    // Validate the complete region first so section assignment is
                    // atomic with respect to invalid or mixed element sets.
                    for (const ID id : *region) {
                        logging::error(id >= 0 && static_cast<std::size_t>(id) < model._data->elements.size() && model._data->elements[static_cast<std::size_t>(id)] != nullptr,
                            "MASS: compiled element ", id, " does not exist");
                        logging::error(model._data->elements[static_cast<std::size_t>(id)]->as<model::PointElement>() != nullptr,
                            "MASS: compiled element set ", *elset, " contains non-MASS element ", id);
                    }

                    // add_section() now receives an assembly-space region and
                    // therefore stores and binds the section directly.
                    model.add_section(std::make_shared<PointMassSection>(region, mass));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
