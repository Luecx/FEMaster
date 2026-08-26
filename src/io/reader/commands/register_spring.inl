/**
 * @file register_spring.inl
 * @brief Registers linear SPRING1 point-element properties.
 *
 * `SPRING, ELSET=...` assigns one scalar ground stiffness to one-node `SPRING1`
 * topology. Degrees of freedom 1--3 map to the translational spring vector of
 * `PointMassSection`; degrees of freedom 4--6 map to its rotational spring
 * vector. The same syntax is available in native FEMaster and Abaqus input.
 *
 * Native FEMaster also retains the legacy `*POINTMASS, NSET=...` variant. It can
 * assign diagonal translational and rotational ground-spring stiffness directly
 * to nodes while optionally combining the spring terms with mass and inertia.
 *
 * The supported subset is intentionally limited to real, constant, linear
 * SPRING1 behavior in the global nodal system. SPRING2/SPRINGA topology,
 * orientations, nonlinear force--displacement data, field/frequency dependence
 * and complex stiffness are not represented by this command.
 *
 * Part-local assignments are created before `Model::compile()` and inherited by
 * Instances. Assembly-level assignments are deferred until compiled ELSETs are
 * available.
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
 * Registers constant linear `*SPRING` properties for one-node SPRING1 elements.
 *
 * Abaqus SPRING1 data consists of the active degree of freedom on the first data
 * line and the scalar stiffness on the second. FEMaster maps DOFs 1--3 to
 * translational stiffness and DOFs 4--6 to rotational stiffness in the existing
 * diagonal point section. Native `*POINTMASS` remains the legacy NSET-based
 * alternative for diagonal ground-spring properties.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Model receiving Part-local or compiled section assignments.
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

        if (!model._data->compiled) {
            // ROOT/PART assignments are resolvable in the active semantic Part.
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
                        // Resolve and validate the complete Part-local point region.
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

                        // Map the Abaqus generalized direction to the diagonal
                        // translational or rotational point-section stiffness.
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

            // Assembly ELSETs do not exist until after compile().
            command.variant(dsl::Variant::make()
                .when(dsl::Condition::parent_is("ASSEMBLY"))
                .segment(dsl::Segment::make()
                    .range(dsl::LineRange{}.min(2).max(2))
                    .pattern(dsl::Pattern::make()
                        .allow_multiline()
                        .one<ID>().name("DOF").desc("Nodal degree of freedom 1..6")
                        .one<Precision>().name("STIFFNESS").desc("Linear spring stiffness")
                    )
                )
            );
            return;
        }

        // ROOT/PART assignments were already copied through compile().
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(2).max(2))
                .pattern(dsl::Pattern::make()
                    .allow_multiline()
                    .one<ID>().name("DOF").desc("Nodal degree of freedom 1..6")
                    .one<Precision>().name("STIFFNESS").desc("Linear spring stiffness")
                )
            )
        );

        // Assembly assignments resolve the compiled point-element region directly.
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
                    // Resolve and validate the complete compiled point region.
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
