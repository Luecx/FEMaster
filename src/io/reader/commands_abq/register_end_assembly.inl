/**
 * @file register_end_assembly.inl
 * @brief Registers the Abaqus ENDASSEMBLY scope terminator.
 *
 * `END ASSEMBLY` closes the syntactic assembly scope and clears the shared
 * parser-local flag used to interpret qualified sets, surfaces and Instances.
 * FEMaster's model-level assembly remains alive because the keyword controls
 * parsing context rather than object lifetime.
 *
 * Registering the terminator in every pass preserves Abaqus nesting even when
 * its callback is inactive and the deck is only being consumed.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_end_assembly(fem::io::dsl::Registry& registry,
                                  std::shared_ptr<bool> assembly_scope) {
    registry.command("ENDASSEMBLY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ASSEMBLY"));
        command.doc("Terminate the Abaqus assembly scope.");
        command.on_enter([assembly_scope](const fem::io::dsl::Keys&) {
            *assembly_scope = false;
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
