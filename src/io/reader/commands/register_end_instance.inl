/**
 * @file register_end_instance.inl
 * @brief Registers the ENDINSTANCE scope terminator.
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
