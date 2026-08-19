/**
 * @file register_assembly.inl
 * @brief Registers the ASSEMBLY scope.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_assembly(fem::io::dsl::Registry& registry,
                              std::shared_ptr<bool> assembly_scope) {
    registry.command("ASSEMBLY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin the assembly scope represented directly by Model.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional("ASSEMBLY").doc("Optional assembly name")
        );
        command.on_enter([assembly_scope](const fem::io::dsl::Keys&) {
            *assembly_scope = true;
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
