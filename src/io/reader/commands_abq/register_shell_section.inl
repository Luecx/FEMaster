/**
 * @file register_shell_section.inl
 * @brief Registers homogeneous Abaqus *SHELL SECTION definitions.
 *
 * The registration assigns one material and one constant thickness to an Abaqus
 * element set and maps the definition to FEMaster's integrated shell section.
 * The supported section uses five Simpson material points through the thickness
 * and may reference a previously defined material orientation.
 *
 * @see model::Model::shell_section
 * @see IntegratedShellSection
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <memory>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers homogeneous Abaqus shell-section assignment.
 *
 * `ELSET` and `MATERIAL` identify the target region and material. `ORIENTATION`
 * optionally selects a material coordinate system. The data line provides the
 * physical shell thickness and an optional integration-point count, which must
 * equal the five Simpson points used by the FEMaster section implementation.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the section definition.
 */
inline void register_shell_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SHELLSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a homogeneous Abaqus shell section to an element set.");

        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET")
                    .required()
                    .doc("Target shell element set")
                .key("MATERIAL")
                    .required()
                    .doc("Material assigned to the homogeneous shell")
                .key("ORIENTATION")
                    .optional()
                    .doc("Optional material orientation")
        );

        command.on_enter([material, elset, orientation](const fem::io::dsl::Keys& keys) {
            *material    = keys.raw("MATERIAL");
            *elset       = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("THICKNESS").desc("Shell thickness")
                    .one<int>().name("INTEGRATIONPOINTS").desc("Through-thickness integration points")
                        .on_missing(5).on_empty(5)
                )
                .bind([&model, material, elset, orientation](fem::Precision thickness, int integration_points) {
                    logging::error(integration_points == 5,
                        "SHELL SECTION supports exactly 5 Simpson integration points");

                    model.shell_section(*elset, *material, thickness, *orientation);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
