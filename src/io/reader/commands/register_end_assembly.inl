/**
 * @file register_end_assembly.inl
 * @brief Registers the ENDASSEMBLY scope terminator.
 *
 * `ENDASSEMBLY` closes the assembly grammar scope and clears the registry-local
 * flag shared with assembly-aware topology commands. The flag determines
 * whether sets, surfaces and other qualified definitions address compiled
 * assembly entities rather than active Part storage.
 *
 * No model object is destroyed when the scope ends because FEMaster represents
 * the assembly directly through the persistent `Model`.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_end_assembly(fem::io::dsl::Registry& registry,
                                  std::shared_ptr<bool> assembly_scope) {
    registry.command("ENDASSEMBLY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ASSEMBLY"));
        command.doc("Terminate the assembly scope.");
        command.on_enter([assembly_scope](const fem::io::dsl::Keys&) {
            *assembly_scope = false;
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
