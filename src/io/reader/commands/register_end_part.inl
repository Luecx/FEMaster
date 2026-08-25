/**
 * @file register_end_part.inl
 * @brief Registers the ENDPART scope terminator.
 *
 * `ENDPART` explicitly terminates an active reusable Part grammar scope. The
 * corresponding model-state transition belongs to the `PART` scope itself and
 * is performed by its exit hook, so explicit and implicit scope exits have the
 * same behavior.
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
        command.doc("Terminate the current part definition.");
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
