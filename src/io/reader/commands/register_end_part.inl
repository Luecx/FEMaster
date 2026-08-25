/**
 * @file register_end_part.inl
 * @brief Registers the ENDPART scope terminator.
 *
 * `ENDPART` is valid only while a `PART` scope is active. The DSL terminator
 * closes that admitted parent directly, which fires the Part exit hook and
 * restores the model's default Part for explicit and implicit scope exits.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_end_part(fem::io::dsl::Registry& registry, model::Model&) {
    registry.command("ENDPART", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("PART"));
        command.closes_parent();
        command.doc("Terminate the current part definition.");
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
