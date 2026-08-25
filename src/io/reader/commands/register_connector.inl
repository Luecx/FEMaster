/**
 * @file register_connector.inl
 * @brief Registers the CONNECTOR input command.
 *
 * `CONNECTOR` relates two compiled nodes through one of the supported
 * `constraint::ConnectorType` masks. Each node is selected through a node set
 * containing exactly one entry, and the referenced coordinate system defines
 * the local directions associated with the connector mask.
 *
 * Registration resolves the named model objects and appends the resulting
 * `constraint::Connector` directly to `ModelData::connectors`. Construction of
 * linear constraint equations remains in `constraint::Connector` and occurs
 * later when `Model::collect_constraints()` assembles the model constraints.
 *
 * @see constraint::Connector
 * @see constraint::ConnectorType
 * @see model::Model::collect_constraints
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <string>

#include "../../../constraints/types/connector.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers two-node connector constraints in the compiled model.
 *
 * The command maps its textual connector type to the corresponding generalized
 * DOF mask, resolves two one-node regions and one coordinate system, and stores
 * the completed connector in model order. It consumes no data lines.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Compiled model providing node regions, coordinate systems and
 *              connector storage.
 */
inline void register_connector(dsl::Registry& registry, model::Model& model) {
    registry.command("CONNECTOR", [&](dsl::Command& command) {
        // Connector endpoints use compiled assembly node sets at root scope
        command.allow_if(dsl::Condition::parent_is("ROOT"));

        // Expose the supported formulations and endpoint requirements in the registry documentation
        command.doc(
            "Relate two one-node sets using a BEAM, HINGE, CYLINDRICAL, "
            "TRANSLATOR, JOIN or JOINRX connector in a named coordinate system."
        );

        // Define the connector formulation, endpoints and local coordinate system
        command.keyword(
            dsl::KeywordSpec::make()
                .key("TYPE").required().doc("Connector formulation")
                .key("NSET1").required().doc("Node set containing the first endpoint")
                .key("NSET2").required().doc("Node set containing the second endpoint")
                .key("COORDINATESYSTEM").required().doc("Coordinate system defining the local connector directions")
        );

        // Resolve and store the connector when its keyword line is entered
        command.on_enter([&model](const dsl::Keys& keys) {
            // Extract the normalized formulation and model references
            const std::string type                   = keys.raw("TYPE");
            const std::string node_set_1_name        = keys.raw("NSET1");
            const std::string node_set_2_name        = keys.raw("NSET2");
            const std::string coordinate_system_name = keys.raw("COORDINATESYSTEM");

            // Map the input formulation to the generalized connector DOF mask
            constraint::ConnectorType connector_type = constraint::ConnectorType::None;
            if (type == "BEAM") {
                connector_type = constraint::ConnectorType::Beam;
            } else if (type == "HINGE") {
                connector_type = constraint::ConnectorType::Hinge;
            } else if (type == "CYLINDRICAL") {
                connector_type = constraint::ConnectorType::Cylindrical;
            } else if (type == "TRANSLATOR") {
                connector_type = constraint::ConnectorType::Translator;
            } else if (type == "JOIN") {
                connector_type = constraint::ConnectorType::Join;
            } else if (type == "JOINRX") {
                connector_type = constraint::ConnectorType::JoinRx;
            }

            // Validate the compiled state, formulation and referenced definitions
            logging::error(model._data->compiled,
                "CONNECTOR: constraints require a compiled model");
            logging::error(connector_type != constraint::ConnectorType::None,
                "CONNECTOR: unsupported type ", type);
            logging::error(model._data->node_sets.has(node_set_1_name),
                "CONNECTOR: node set ", node_set_1_name, " does not exist");
            logging::error(model._data->node_sets.has(node_set_2_name),
                "CONNECTOR: node set ", node_set_2_name, " does not exist");
            logging::error(model._data->coordinate_systems.has(coordinate_system_name),
                "CONNECTOR: coordinate system ", coordinate_system_name, " does not exist");

            // Resolve the validated references once for cardinality checks and construction
            const auto node_set_1        = model._data->node_sets.get(node_set_1_name);
            const auto node_set_2        = model._data->node_sets.get(node_set_2_name);
            const auto coordinate_system = model._data->coordinate_systems.get(coordinate_system_name);

            logging::error(node_set_1 != nullptr && node_set_1->size() == 1,
                "CONNECTOR: node set ", node_set_1_name, " must contain exactly one node");
            logging::error(node_set_2 != nullptr && node_set_2->size() == 1,
                "CONNECTOR: node set ", node_set_2_name, " must contain exactly one node");

            // Store the connector; equation construction is deferred to constraint collection
            model._data->connectors.emplace_back(
                node_set_1->first(),
                node_set_2->first(),
                coordinate_system,
                connector_type
            );
        });

        // CONNECTOR is fully defined by its keyword arguments and has no data payload
        command.variant(dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
