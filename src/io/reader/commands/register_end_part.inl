/**
 * @file register_end_part.inl
 * @brief Registers the ENDPART scope terminator.
 *
 * `ENDPART` closes the active reusable Part definition and restores the model's
 * default part as the target for subsequent unqualified topology commands. The
 * transition is permitted only before model compilation because semantic Part
 * topology becomes immutable at the compile boundary.
 *
 * Scope admission is enforced by the DSL registry; this callback performs the
 * corresponding model-state transition.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_end_part(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ENDPART", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("PART"));
        command.doc("Terminate the current part definition.");
        command.on_enter([&model](const fem::io::dsl::Keys&) {
            logging::error(model._data != nullptr && !model._data->compiled,
                "ENDPART: part scope cannot change after compile()");
            model._data->parts.activate(model::Model::DEFAULT_PART_NAME);
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
