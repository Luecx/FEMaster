/**
 * @file register_mpc.cpp
 * @brief Registers simple Abaqus-style multi-point constraints.
 *
 * The shared MPC command accepts the two-node BEAM, TIE and PIN relations used
 * by both native and Abaqus input modes. Each data line is resolved after model
 * compilation and stored as one constraint::Mpc definition.
 *
 * Equation generation remains inside constraint::Mpc and is deferred until the
 * active structural DOF map is built for a load case.
 *
 * @see constraint::Mpc
 * @see model::Model::collect_constraints
 *
 * @author Finn Eggers
 * @date 01.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <string>

#include "../../../constraints/types/mpc.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * Registers the MPC command for compiled root and assembly scopes.
 *
 * Every line uses Abaqus ordering
 *
 *     TYPE, dependent node, independent node
 *
 * and creates one persistent Mpc object. Only the simple two-node BEAM, TIE and
 * PIN types are accepted; more specialized MPC formulations require additional
 * geometric definitions and are intentionally outside this command.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Compiled model receiving the MPC definitions.
 */
void register_mpc(dsl::Registry& registry, model::Model& model) {
    registry.command("MPC", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Define BEAM, TIE or PIN multi-point constraints between two nodes.");

        command.on_enter([&model](const dsl::Keys&) {
            logging::error(model._data->compiled,
                "MPC: constraints require a compiled model");
        });

        command.data(
            dsl::Pattern::make()
                .one<std::string>("TYPE")
                .one<std::string>("NODE1")
                .one<std::string>("NODE2"),
            [&model](const std::string& type,
                           const std::string& node_1,
                           const std::string& node_2) {
                constraint::Mpc::Type mpc_type = constraint::Mpc::Type::Beam;
                if (type == "BEAM") {
                    mpc_type = constraint::Mpc::Type::Beam;
                } else if (type == "TIE") {
                    mpc_type = constraint::Mpc::Type::Tie;
                } else if (type == "PIN") {
                    mpc_type = constraint::Mpc::Type::Pin;
                } else {
                    logging::error(false,
                        "MPC: unsupported type ", type, "; supported types are BEAM, TIE and PIN");
                }

                model._data->mpcs.emplace_back(
                    model.compiled_node_id(node_1),
                    model.compiled_node_id(node_2),
                    mpc_type
                );
            }
        );
    });
}

} // namespace fem::io::reader::commands
