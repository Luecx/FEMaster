/**
 * @file register_end_assembly.inl
 * @brief Registers the ENDASSEMBLY scope terminator.
 *
 * `ENDASSEMBLY` closes the assembly grammar scope. FEMaster does not destroy or
 * modify a separate assembly object because the assembled domain is represented
 * directly by the persistent `Model`.
 *
 * Scope ownership remains entirely inside the DSL engine; no parser-local flag
 * needs to be cleared explicitly.
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
 * @brief Registers the `ENDASSEMBLY` grammar-scope terminator.
 */
inline void register_end_assembly(dsl::Registry& registry) {
    registry.command("ENDASSEMBLY", [&](dsl::Command& command) {
        // ENDASSEMBLY is only valid for an active assembly scope
        command.allow_if(dsl::Condition::parent_is("ASSEMBLY"));
        command.closes_parent();

        // The command only closes the grammar scope
        command.variant(dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
