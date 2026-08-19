/**
 * @file register_connector.inl
 * @brief Registers connector constraints.
 *
 * The parser resolves all referenced model objects and constructs the concrete
 * Connector itself. Model deliberately has no generic constraint dispatcher;
 * the completed object is appended directly to ModelData.
 *
 * The command binds two node regions through a selected connector formulation
 * and coordinate system. The resulting constraint contributes its kinematic
 * relations when an active load case collects model constraints.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../../constraints/types/connector.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_connector(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONNECTOR", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").required()
                .key("NSET1").required()
                .key("NSET2").required()
                .key("COORDINATESYSTEM").required()
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string type = keys.raw("TYPE");
            const std::string set1 = keys.raw("NSET1");
            const std::string set2 = keys.raw("NSET2");
            const std::string csys = keys.raw("COORDINATESYSTEM");

            logging::error(model._data->compiled,
                "CONNECTOR: constraints require a compiled model");
            logging::error(model._data->node_sets.has(set1),
                "CONNECTOR: node set ", set1, " does not exist");
            logging::error(model._data->node_sets.get(set1) && model._data->node_sets.get(set1)->size() == 1,
                "CONNECTOR: node set ", set1, " must contain exactly one node");
            logging::error(model._data->node_sets.has(set2),
                "CONNECTOR: node set ", set2, " does not exist");
            logging::error(model._data->node_sets.get(set2) && model._data->node_sets.get(set2)->size() == 1,
                "CONNECTOR: node set ", set2, " must contain exactly one node");
            logging::error(model._data->coordinate_systems.has(csys),
                "CONNECTOR: coordinate system ", csys, " does not exist");

            constraint::ConnectorType connector_type = constraint::ConnectorType::None;
            if (type == "BEAM") connector_type = constraint::ConnectorType::Beam;
            else if (type == "HINGE") connector_type = constraint::ConnectorType::Hinge;
            else if (type == "CYLINDRICAL") connector_type = constraint::ConnectorType::Cylindrical;
            else if (type == "TRANSLATOR") connector_type = constraint::ConnectorType::Translator;
            else if (type == "JOIN") connector_type = constraint::ConnectorType::Join;
            else if (type == "JOINRX") connector_type = constraint::ConnectorType::JoinRx;

            logging::error(connector_type != constraint::ConnectorType::None,
                "CONNECTOR: unsupported type ", type);

            model._data->connectors.emplace_back(
                model._data->node_sets.get(set1)->first(),
                model._data->node_sets.get(set2)->first(),
                model._data->coordinate_systems.get(csys),
                connector_type
            );
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
