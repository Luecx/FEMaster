/**
 * @file register_shell_section.inl
 * @brief Registers homogeneous Abaqus *SHELL SECTION definitions.
 *
 * The supported Abaqus form assigns one material and one constant thickness to
 * an element set. It maps directly onto `model::Model::shell_section()` and the
 * existing `IntegratedShellSection`, which evaluates the material at five
 * Simpson points through the physical thickness.
 *
 * Abaqus also defaults to five Simpson section points for a homogeneous shell.
 * An explicitly supplied integration-point count is therefore accepted only
 * when it is five. Composite sections, variable thickness definitions, offsets,
 * additional section density and alternative integration rules are intentionally
 * left unsupported.
 *
 * @see model::Model::shell_section
 * @see IntegratedShellSection
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers homogeneous Abaqus shell-section assignment.
 *
 * The keyword requires `ELSET` and `MATERIAL` and optionally references an
 * already defined material orientation. The single data line contains the shell
 * thickness and an optional through-thickness integration-point count. FEMaster
 * uses five Simpson material points internally, matching the Abaqus default, so
 * any explicitly supplied value other than five is rejected rather than silently
 * changing the requested section integration.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the section definition.
 */
inline void register_shell_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SHELLSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a homogeneous Abaqus shell section to an element set.");

        // Preserve the keyword data until the following section data line is
        // parsed and the complete FEMaster section can be constructed
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

        // Read the physical thickness and retain the common Abaqus/FEMaster
        // default of five Simpson points through the shell section
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("THICKNESS").desc("Shell thickness")
                    .one<int>().name("INTEGRATIONPOINTS").desc("Through-thickness integration points")
                        .on_missing(5).on_empty(5)
                )
                .bind([&model, material, elset, orientation](fem::Precision thickness, int integration_points) {
                    if (integration_points != 5) {
                        throw std::runtime_error("SHELL SECTION currently supports exactly 5 Simpson integration points");
                    }

                    model.shell_section(*elset, *material, thickness, *orientation);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
