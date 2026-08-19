/**
 * @file register_contact.inl
 * @brief Registers frictionless surface contact.
 *
 * Contact syntax is fully resolved by the parser. Once the master/slave surface
 * sets and penalty parameters are validated, the concrete Contact object is
 * appended directly to ModelData. Model itself remains free of parser-oriented
 * constraint factories and type dispatch.
 *
 * The stored definition includes search distance, penalty stiffness, clearance
 * and optional master-normal reversal. Contact search, active-point updates and
 * force/tangent assembly remain inside the concrete contact implementation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../../constraints/types/contact.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_contact(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONTACT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER").required()
                .key("SLAVE").required()
                .key("PENALTY").required()
                .key("CLEARANCE").optional("0")
                .key("FLIP").optional("NO").allowed({"NO", "YES"})
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string master = keys.raw("MASTER");
            const std::string slave  = keys.raw("SLAVE");

            logging::error(model._data->compiled,
                "CONTACT: constraints require a compiled model");
            logging::error(model._data->surface_sets.has(master),
                "CONTACT: master surface set ", master, " does not exist");
            logging::error(model._data->surface_sets.get(master) && model._data->surface_sets.get(master)->size() > 0,
                "CONTACT: master surface set ", master, " is empty");
            logging::error(model._data->surface_sets.has(slave),
                "CONTACT: slave surface set ", slave, " does not exist");
            logging::error(model._data->surface_sets.get(slave) && model._data->surface_sets.get(slave)->size() > 0,
                "CONTACT: slave surface set ", slave, " is empty");

            model._data->contacts.emplace_back(
                model._data->surface_sets.get(master),
                model._data->surface_sets.get(slave),
                keys.get<Precision>("PENALTY"),
                keys.get<Precision>("CLEARANCE"),
                keys.get<bool>("FLIP")
            );
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
