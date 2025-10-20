// register_connector.inl — registers *CONNECTOR
#pragma once
/**
 * @file register_connector.inl
 * @brief Register *CONNECTOR command.
 *
 * Defines a connector constraint between two node sets using a named coordinate system.
 */

#include <string>

#include "../../constraints/connector.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"
#include "../../core/logging.h" // logging::error(...)

namespace fem::input_decks::commands {

inline void register_connector(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("CONNECTOR", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));

        command.doc(
            "Define a connector between two node sets in a given coordinate system. "
            "Both node sets must exist and be non-empty. Their composition should fit the "
            "connector’s semantics. The coordinate system defines the local axes."
        );

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE").required().doc("Connector type identifier.")
                .key("NSET1").required().doc("First node set.")
                .key("NSET2").required().doc("Second node set.")
                .key("COORDINATESYSTEM").required().doc("Name of a coordinate system.")
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            const std::string type_raw = keys.raw("TYPE");   // already uppercase per DSL
            const std::string nset1    = keys.raw("NSET1");
            const std::string nset2    = keys.raw("NSET2");
            const std::string coord    = keys.raw("COORDINATESYSTEM");

            // Existence
            if (!model._data->node_sets.has(nset1))
                logging::error(false, "CONNECTOR: node set '" + nset1 + "' does not exist.");
            if (!model._data->node_sets.has(nset2))
                logging::error(false, "CONNECTOR: node set '" + nset2 + "' does not exist.");

            // Emptiness (Region ptr -> use ->size())
            if (model._data->node_sets.has(nset1) && model._data->node_sets.get(nset1)->size() == 0)
                logging::error(false, "CONNECTOR: node set '" + nset1 + "' is empty.");
            if (model._data->node_sets.has(nset2) && model._data->node_sets.get(nset2)->size() == 0)
                logging::error(false, "CONNECTOR: node set '" + nset2 + "' is empty.");

            // TYPE mapping (inputs already uppercase)
            constraint::ConnectorType ctype = constraint::ConnectorType::None;
            if      (type_raw == "BEAM")        ctype = constraint::ConnectorType::Beam;
            else if (type_raw == "HINGE")       ctype = constraint::ConnectorType::Hinge;
            else if (type_raw == "CYLINDRICAL") ctype = constraint::ConnectorType::Cylindrical;
            else if (type_raw == "TRANSLATOR")  ctype = constraint::ConnectorType::Translator;
            else if (type_raw == "JOIN")        ctype = constraint::ConnectorType::Join;
            else if (type_raw == "JOINRX")      ctype = constraint::ConnectorType::JoinRx;
            else
                logging::error(false,
                               "CONNECTOR: unknown TYPE='" + type_raw +
                               "'. Supported: BEAM, HINGE, CYLINDRICAL, TRANSLATOR, JOIN, JOINRX.");

            // Delegate (coord-system checks handled inside model)
            model.add_connector(nset1, nset2, coord, ctype);
        });

        command.variant(fem::dsl::Variant::make()); // no data lines
    });
}

} // namespace fem::input_decks::commands
