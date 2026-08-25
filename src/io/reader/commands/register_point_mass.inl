/**
 * @file register_point_mass.inl
 * @brief Registers FEMaster point-mass sections on POINT element sets.
 *
 * `POINTMASS` is a section/property command rather than a nodal feature. The
 * target ELSET must contain only one-node `PointElement` topology. The section
 * supplies translational mass, rotary inertia and diagonal translational and
 * rotational stiffness/damping against ground.
 *
 * Scope selection is declarative through `parent_is(...)` on the variant; the
 * bind contains no parser-stage or scope branching.
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../model/element/point.h"
#include "../../../model/model.h"
#include "../../../section/section_point_mass.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers a concentrated point section on an element set.
 *
 * Data layout remains compact and backward-compatible in its first ten values:
 * scalar mass, three rotary inertias, three translational ground stiffnesses and
 * three rotational ground stiffnesses. Two optional three-vectors append
 * translational and rotational viscous ground damping.
 *
 * @param registry Parser registry receiving the command.
 * @param model Model receiving the part-local section assignment.
 */
inline void register_point_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("POINTMASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));
        command.doc("Assign concentrated mass, inertia and ground impedance to POINT elements.");

        auto elset = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").required().doc("Target POINT element set")
        );

        command.on_enter([elset](const fem::io::dsl::Keys& keys) {
            *elset = keys.raw("ELSET");
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::Precision>().name("MASS").desc("Point mass")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("INERTIA").desc("Rotary inertia Ix,Iy,Iz")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("SPRING").desc("Translational ground stiffness")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("ROTSPRING").desc("Rotational ground stiffness")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("DAMPING").desc("Translational ground damping")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("ROTDAMPING").desc("Rotational ground damping")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, elset](
                          fem::Precision mass,
                          const std::array<fem::Precision, 3>& inertia_data,
                          const std::array<fem::Precision, 3>& spring_data,
                          const std::array<fem::Precision, 3>& rotary_spring_data,
                          const std::array<fem::Precision, 3>& damping_data,
                          const std::array<fem::Precision, 3>& rotary_damping_data) {
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "POINTMASS: no active part is available");
                    logging::error(part->elem_sets.has(*elset),
                        "POINTMASS: element set ", *elset, " does not exist in part ", part->name);

                    const auto region = part->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "POINTMASS: element set ", *elset, " is empty in part ", part->name);

                    for (const ID id : *region) {
                        const auto it = part->elements.find(id);
                        logging::error(it != part->elements.end() && it->second != nullptr,
                            "POINTMASS: element ", id, " does not exist in part ", part->name);
                        logging::error(it->second->as<model::PointElement>() != nullptr,
                            "POINTMASS: element set ", *elset, " contains non-POINT element ", id);
                    }

                    const fem::Vec3 inertia{inertia_data[0], inertia_data[1], inertia_data[2]};
                    const fem::Vec3 spring{spring_data[0], spring_data[1], spring_data[2]};
                    const fem::Vec3 rotary_spring{rotary_spring_data[0], rotary_spring_data[1], rotary_spring_data[2]};
                    const fem::Vec3 damping{damping_data[0], damping_data[1], damping_data[2]};
                    const fem::Vec3 rotary_damping{rotary_damping_data[0], rotary_damping_data[1], rotary_damping_data[2]};

                    model.add_section(std::make_shared<PointMassSection>(
                        region, mass, inertia, spring, rotary_spring, damping, rotary_damping));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
