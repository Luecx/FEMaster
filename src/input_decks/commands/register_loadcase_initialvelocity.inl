#pragma once
/**
 * @file register_loadcase_initialvelocity.inl
 * @brief Register *INITIALVELOCITY inside *LOADCASE (Transient): reference a node field.
 */

#include <string>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../core/logging.h"
#include "../../data/field.h"

#include "../../loadcase/linear_transient.h"

namespace fem::input_decks::commands {

inline void register_loadcase_initialvelocity(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("INITIALVELOCITY", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set initial velocity for transient analysis from a node FIELD. Usage: *INITIALVELOCITY, FIELD=NAME");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("FIELD").required().doc("Name of a node field with 6 components (ux,uy,uz,rx,ry,rz)")
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            logging::error(base != nullptr, "INITIALVELOCITY must appear inside *LOADCASE.");

            const std::string field = keys.raw("FIELD");

            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                // Resolve field now and attach pointer directly
                auto& mdl = parser.model();
                auto fptr = mdl._data->get_field(field);
                logging::error(fptr != nullptr, "Initial velocity field ", field, " does not exist");
                logging::error(fptr->domain == fem::model::FieldDomain::NODE, "Initial velocity field ", field, " must be a node field");
                logging::error(fptr->components == 6, "Initial velocity field ", field, " must have exactly 6 components");
                lc->initial_velocity = fptr;
                return;
            }

            logging::error(false, "INITIALVELOCITY not supported for loadcase type " + parser.active_loadcase_type());
        });

        // No data lines for this command
        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
