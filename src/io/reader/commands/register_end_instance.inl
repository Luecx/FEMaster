/**
 * @file register_end_instance.inl
 * @brief Registers the ENDINSTANCE scope terminator.
 *
 * `ENDINSTANCE` closes the syntactic scope opened by an `INSTANCE` command.
 * Instance placement is already captured when the opening command executes, so
 * no additional model mutation is required at the terminator.
 *
 * Keeping the terminator in the complete grammar is nevertheless essential for
 * correct nested-scope consumption during every dependency-ordered parser pass.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_end_instance(fem::io::dsl::Registry& registry) {
    registry.command("ENDINSTANCE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("INSTANCE"));
        command.doc("Terminate the current instance scope.");
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
