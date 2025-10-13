// register_connector.inl — registers *CONNECTOR

#include <stdexcept>
#include <string>

#include "../../constraints/connector.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_connector(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("CONNECTOR", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Define a connector between two node sets.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE").required().doc("Connector type")
                .key("NSET1").required().doc("First node set")
                .key("NSET2").required().doc("Second node set")
                .key("COORDINATESYSTEM").required().doc("Coordinate system name")
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            const std::string& type = keys.raw("TYPE");
            const std::string& nset1 = keys.raw("NSET1");
            const std::string& nset2 = keys.raw("NSET2");
            const std::string& coord = keys.raw("COORDINATESYSTEM");

            constraint::ConnectorType ctype = constraint::ConnectorType::None;
            if (type == "BEAM") ctype = constraint::ConnectorType::Beam;
            else if (type == "HINGE") ctype = constraint::ConnectorType::Hinge;
            else if (type == "CYLINDRICAL") ctype = constraint::ConnectorType::Cylindrical;
            else if (type == "TRANSLATOR") ctype = constraint::ConnectorType::Translator;
            else if (type == "JOIN") ctype = constraint::ConnectorType::Join;
            else if (type == "JOINRX") ctype = constraint::ConnectorType::JoinRx;
            else throw std::runtime_error("Unknown connector TYPE=" + type);

            model.add_connector(nset1, nset2, coord, ctype);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

