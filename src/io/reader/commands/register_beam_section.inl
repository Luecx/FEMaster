/**
 * @file register_beam_section.inl
 * @brief Registers beam section assignments.
 *
 * `BEAMSECTION` combines a material, beam profile and target element region into
 * a `BeamSection`. An optional orientation vector defines the initial local
 * cross-section direction used to construct the beam frame.
 *
 * The section is attached to semantic Part topology before compilation. Element
 * kinematics and constitutive integration consume the stored profile, material
 * and orientation later during matrix and result assembly.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>

#include "../../../model/model.h"
#include "../../../section/section_beam.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_beam_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("BEAMSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));

        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto profile     = std::make_shared<std::string>();
        auto orientation = std::make_shared<Vec3>(Vec3::Zero());

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required()
                .key("ELSET").required()
                .key("PROFILE").required()
        );
        command.on_enter([material, elset, profile, orientation](const fem::io::dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset    = keys.raw("ELSET");
            *profile  = keys.raw("PROFILE");
            orientation->setZero();
        });
        command.on_exit([&model, material, elset, profile, orientation](const fem::io::dsl::Keys&) {
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "BEAMSECTION: no active part is available");
            logging::error(part->elem_sets.has(*elset),
                "BEAMSECTION: element set ", *elset, " is not defined in part ", part->name);
            logging::error(model._data->materials.has(*material),
                "BEAMSECTION: material ", *material, " is not defined");
            logging::error(model._data->profiles.has(*profile),
                "BEAMSECTION: profile ", *profile, " is not defined");
            logging::error(orientation->norm() > Precision(0),
                "BEAMSECTION: orientation vector must be non-zero");

            auto section = std::make_shared<BeamSection>();
            section->material_  = model._data->materials.get(*material);
            section->region_    = part->elem_sets.get(*elset);
            section->profile_   = model._data->profiles.get(*profile);
            section->direction_ = *orientation;
            model.add_section(std::move(section));
        });
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<Precision, 3>().name("N1")
                )
                .bind([orientation](const std::array<Precision, 3>& values) {
                    *orientation = Vec3{values[0], values[1], values[2]};
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
