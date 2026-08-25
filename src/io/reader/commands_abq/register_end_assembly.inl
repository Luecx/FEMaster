/**
 * @file register_end_assembly.inl
 * @brief Registers the Abaqus ENDASSEMBLY scope terminator.
 *
 * `ENDASSEMBLY` closes the Abaqus assembly grammar scope. FEMaster's model-level
 * assembly remains alive because the keyword controls parsing context rather
 * than object lifetime.
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
 * @brief Registers the Abaqus `ENDASSEMBLY` grammar-scope terminator.
 */
inline void register_end_assembly(dsl::Registry& registry) {
    registry.command("ENDASSEMBLY", [&](dsl::Command& command) {
        // ENDASSEMBLY is only valid for an active assembly scope
        command.allow_if(dsl::Condition::parent_is("ASSEMBLY"));

        // The command only closes the grammar scope
        command.variant(dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
