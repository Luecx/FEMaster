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
 */
inline void register_assembly(dsl::Registry& registry) {
    registry.command("ASSEMBLY", [&](dsl::Command& command) {
        // Restrict ASSEMBLY to the root scope
        command.allow_if(dsl::Condition::parent_is("ROOT"));

        // Define the optional assembly name accepted by the input format
        command.keyword(
            dsl::KeywordSpec::make()
                .key("NAME").optional("ASSEMBLY").doc("Optional assembly name")
        );

        // ASSEMBLY only introduces a grammar scope
        command.variant(dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
