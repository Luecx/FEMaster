/**
 * @file register_assembly.inl
 * @brief Registers the Abaqus ASSEMBLY scope.
 *
 * FEMaster has one model-level assembly, so the Abaqus `ASSEMBLY` keyword only
 * opens the corresponding grammar scope. Parent-dependent set and surface
 * behavior is selected directly through DSL variant conditions.
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers the root-level Abaqus `ASSEMBLY` grammar scope.
 */
inline void register_assembly(dsl::Registry& registry) {
    registry.command("ASSEMBLY", [&](dsl::Command& command) {
        // Restrict ASSEMBLY to the root scope
        command.allow_if(dsl::Condition::parent_is("ROOT"));

        // Define the optional Abaqus assembly name
        command.keyword(
            dsl::KeywordSpec::make()
                .key("NAME").optional("ASSEMBLY").doc("Assembly name")
        );

        // ASSEMBLY only introduces a grammar scope
        command.variant(dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
