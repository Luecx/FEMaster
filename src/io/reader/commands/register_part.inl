/**
 * @file register_part.inl
 * @brief Registers the PART scope.
 *
 * A `PART` command creates and activates a named reusable semantic topology
 * container. Subsequent nodes, elements, local regions and section assignments
 * are stored against that Part until `ENDPART` restores the default context.
 *
 * Parts retain sparse local identifiers and remain independent of assembly
 * placement. Rigid transforms and dense global expansion are handled later by
 * Instance construction and `Model::compile()`.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_part(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PART", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin a reusable finite-element part definition.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Unique part name")
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model.add_part(keys.raw("NAME"));
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
