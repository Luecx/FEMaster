/**
 * @file register_assembly.inl
 * @brief Registers the Abaqus ASSEMBLY scope.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_assembly(fem::io::dsl::Registry& registry,
                              std::shared_ptr<bool> assembly_scope) {
    registry.command("ASSEMBLY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin the Abaqus assembly scope.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional("ASSEMBLY").doc("Assembly name")
        );
        command.on_enter([assembly_scope](const fem::io::dsl::Keys&) {
            *assembly_scope = true;
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
