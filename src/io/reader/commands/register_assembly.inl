/**
 * @file register_assembly.inl
 * @brief Registers the ASSEMBLY input scope.
 *
 * FEMaster represents the assembly directly through `Model` and therefore does
 * not create a separate assembly object. `ASSEMBLY` only introduces a grammar
 * scope under which assembly-specific command variants may be selected.
 *
 * Scope ownership remains entirely inside the DSL engine; no parallel parser
 * flag is required.
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers the root-level `ASSEMBLY` grammar scope.
 *
 * The command accepts the optional assembly name used by the input format but
 * does not create or activate a separate model object. Its admitted DSL scope
 * provides the parent context required by assembly-level child commands.
 *
 * @param registry Parser registry receiving the command definition.
 */
inline void register_assembly(dsl::Registry& registry) {
    registry.command("ASSEMBLY", [&](dsl::Command& command) {
        // Restrict ASSEMBLY to the root scope
        command.allow_if(dsl::Condition::parent_is("ROOT"));

        // Describe the syntactic scope in the generated command documentation
        command.doc("Begin the model-level assembly scope.");

        // Accept the optional input name without creating a separate assembly object
        command.keyword(
            dsl::KeywordSpec::make()
                .key("NAME").optional("ASSEMBLY").doc("Optional input-deck assembly name")
        );

        // Declare an empty payload because ASSEMBLY only introduces a grammar scope
        command.variant(dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
