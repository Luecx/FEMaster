/**
 * @file register_end_instance.inl
 * @brief Registers the Abaqus ENDINSTANCE scope terminator.
 *
 * `END INSTANCE` closes the scope opened by an Abaqus `INSTANCE` command. The
 * referenced Part and rigid placement are already stored by the opening
 * callback, so the terminator performs no additional geometry mutation.
 *
 * Its grammar entry remains necessary in every parser pass to retain valid
 * parent-child scope state while commands are replayed consume-only.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_end_instance(fem::io::dsl::Registry& registry) {
    registry.command("ENDINSTANCE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("INSTANCE"));
        command.doc("Terminate the current Abaqus instance scope.");
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
