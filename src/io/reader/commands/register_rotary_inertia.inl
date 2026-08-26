/**
 * @file register_rotary_inertia.inl
 * @brief Registers concentrated rotary-inertia properties for point elements.
 *
 * `ROTARY INERTIA, ELSET=...` assigns the diagonal rotary inertia of one-node
 * `ROTARYI` topology. Both the native FEMaster and supported Abaqus readers use
 * the same command and represent the property through `PointMassSection` on the
 * existing `model::PointElement` implementation.
 *
 * Native FEMaster also retains the legacy `*POINTMASS, NSET=...` variant. It can
 * assign the same diagonal rotary inertia directly to nodes while optionally
 * combining it with translational mass and diagonal ground-spring stiffness.
 *
 * Abaqus permits a full symmetric inertia tensor. FEMaster's current point
 * section stores only the three diagonal moments, so the supported subset accepts
 * all six Abaqus data fields but requires the products of inertia I12, I13 and
 * I23 to be zero. Orientation-dependent inertia is not part of this mapping.
 *
 * As with `MASS`, Part-local assignments are created before `Model::compile()`
 * and copied through Instances. Assembly-level assignments are deferred until
 * the compiled ELSET namespace exists.
 *
 * @see model::PointElement
 * @see PointMassSection
 * @see commands::register_mass
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

#include <array>
#include <memory>
#include <string>

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * Registers diagonal `*ROTARY INERTIA` properties on point-element ELSETs.
 *
 * The first three data values are mapped to the rotational entries of
 * `PointMassSection`. The optional off-diagonal Abaqus tensor entries are parsed
 * for syntax compatibility and rejected unless they are zero because the
 * current point element assembles only a diagonal rotational mass matrix.
 * Native `*POINTMASS` remains the legacy NSET-based alternative.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Model receiving Part-local or compiled section assignments.
 */
inline void register_rotary_inertia(dsl::Registry& registry, model::Model& model) {
    registry.command("ROTARYINERTIA", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Assign diagonal rotary inertia to ROTARYI point elements; native legacy alternative: POINTMASS.");

        auto elset = std::make_shared<std::string>();

        command.keyword(
            dsl::KeywordSpec::make()
                .key("ELSET").required().doc("Element set containing ROTARYI point elements")
        );

        command.on_enter([elset](const dsl::Keys& keys) {
            *elset = keys.raw("ELSET");
        });

        if (!model._data->compiled) {
            // ROOT/PART assignments are resolvable in the active semantic Part.
            command.variant(dsl::Variant::make()
                .when(dsl::Condition::parent_is({"ROOT", "PART"}))
                .segment(dsl::Segment::make()
                    .range(dsl::LineRange{}.min(1).max(1))
                    .pattern(dsl::Pattern::make()
                        .fixed<Precision, 6>().name("INERTIA").desc("I11,I22,I33,I12,I13,I23")
                            .on_missing(Precision(0)).on_empty(Precision(0))
                    )
                    .bind([&model, elset](const std::array<Precision, 6>& values) {
                        // Resolve and validate the complete Part-local point region.
                        const auto part = model._data->parts.get();
                        logging::error(part != nullptr,
                            "ROTARY INERTIA: no active part is available");
                        logging::error(part->elem_sets.has(*elset),
                            "ROTARY INERTIA: element set ", *elset,
                            " does not exist in part ", part->name);

                        const auto region = part->elem_sets.get(*elset);
                        logging::error(region != nullptr && region->size() > 0,
                            "ROTARY INERTIA: element set ", *elset,
                            " is empty in part ", part->name);

                        for (const ID id : *region) {
                            const auto it = part->elements.find(id);
                            logging::error(it != part->elements.end() && it->second != nullptr,
                                "ROTARY INERTIA: element ", id,
                                " does not exist in part ", part->name);
                            logging::error(it->second->as<model::PointElement>() != nullptr,
                                "ROTARY INERTIA: element set ", *elset,
                                " contains non-point element ", id);
                        }

                        // Only diagonal rotary inertia is represented by PointMassSection.
                        logging::error(values[3] == Precision(0)
                                    && values[4] == Precision(0)
                                    && values[5] == Precision(0),
                            "ROTARY INERTIA: products of inertia I12, I13 and I23 are not supported");

                        const Vec3 inertia{values[0], values[1], values[2]};
                        model.add_section(std::make_shared<PointMassSection>(
                            region, Precision(0), inertia));
                    })
                )
            );

            // Assembly ELSETs do not exist until after compile().
            command.variant(dsl::Variant::make()
                .when(dsl::Condition::parent_is("ASSEMBLY"))
                .segment(dsl::Segment::make()
                    .range(dsl::LineRange{}.min(1).max(1))
                    .pattern(dsl::Pattern::make()
                        .fixed<Precision, 6>().name("INERTIA").desc("I11,I22,I33,I12,I13,I23")
                            .on_missing(Precision(0)).on_empty(Precision(0))
                    )
                )
            );
            return;
        }

        // ROOT/PART assignments were already copied through compile().
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1).max(1))
                .pattern(dsl::Pattern::make()
                    .fixed<Precision, 6>().name("INERTIA").desc("I11,I22,I33,I12,I13,I23")
                        .on_missing(Precision(0)).on_empty(Precision(0))
                )
            )
        );

        // Assembly assignments resolve the compiled point-element region directly.
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::parent_is("ASSEMBLY"))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1).max(1))
                .pattern(dsl::Pattern::make()
                    .fixed<Precision, 6>().name("INERTIA").desc("I11,I22,I33,I12,I13,I23")
                        .on_missing(Precision(0)).on_empty(Precision(0))
                )
                .bind([&model, elset](const std::array<Precision, 6>& values) {
                    // Resolve and validate the complete compiled point region.
                    logging::error(model._data->elem_sets.has(*elset),
                        "ROTARY INERTIA: compiled element set ", *elset, " does not exist");

                    const auto region = model._data->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "ROTARY INERTIA: compiled element set ", *elset, " is empty");

                    for (const ID id : *region) {
                        logging::error(id >= 0
                                    && static_cast<std::size_t>(id) < model._data->elements.size()
                                    && model._data->elements[static_cast<std::size_t>(id)] != nullptr,
                            "ROTARY INERTIA: compiled element ", id, " does not exist");
                        logging::error(model._data->elements[static_cast<std::size_t>(id)]->as<model::PointElement>() != nullptr,
                            "ROTARY INERTIA: compiled element set ", *elset,
                            " contains non-point element ", id);
                    }

                    logging::error(values[3] == Precision(0)
                                && values[4] == Precision(0)
                                && values[5] == Precision(0),
                        "ROTARY INERTIA: products of inertia I12, I13 and I23 are not supported");

                    const Vec3 inertia{values[0], values[1], values[2]};
                    model.add_section(std::make_shared<PointMassSection>(
                        region, Precision(0), inertia));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
