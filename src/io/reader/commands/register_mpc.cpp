/**
 * @file register_mpc.cpp
 * @brief Registers simple Abaqus-style multi-point constraints.
 *
 * The shared MPC command accepts the two-node BEAM, TIE and PIN relations used
 * by both native and Abaqus input modes. Node sets are expanded after model
 * compilation into individual constraint::Mpc definitions.
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
 * and creates one Mpc object per resolved node pair. A first-operand NSET
 * expands against a single master node or pairs with an equally sized second
 * NSET in compiled set order. Only the simple two-node BEAM, TIE and
 * PIN types are accepted; more specialized MPC formulations require additional
 * geometric definitions and are intentionally outside this command.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Compiled model receiving the MPC definitions.
 */
void register_mpc(dsl::Registry& registry, model::Model& model) {
    registry.command("MPC", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Define BEAM, TIE or PIN multi-point constraints between two nodes.");

        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TYPE")
                    .one<std::string>().name("NODE1")
                    .one<std::string>().name("NODE2")
                )
                .bind([&model](const std::string& type,
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

                    // Part-scoped MPCs run before compile() and stay part-local
                    const bool local     = !model._data->compiled;
                    const auto part      = model._data->parts.get();
                    auto&      node_sets = local ? part->node_sets : model._data->node_sets;
                    auto&      mpcs      = local ? part->mpcs      : model._data->mpcs;

                    const bool first_is_set  = node_sets.has(node_1);
                    const bool second_is_set = node_sets.has(node_2);

                    logging::error(!second_is_set || first_is_set,
                        "MPC: second node set requires a first node set");

                    if (!first_is_set) {
                        mpcs.emplace_back(
                            model.resolve_node_id(node_1),
                            model.resolve_node_id(node_2),
                            mpc_type
                        );
                        return;
                    }

                    const auto first = node_sets.get(node_1);
                    logging::error(first != nullptr && first->size() > 0,
                        "MPC: first node set ", node_1, " is empty");

                    if (second_is_set) {
                        const auto second = node_sets.get(node_2);
                        logging::error(second != nullptr && second->size() == first->size(),
                            "MPC: node set ", node_2, " has incompatible size");

                        for (std::size_t i = 0; i < first->size(); ++i) {
                            mpcs.emplace_back(first->at(i), second->at(i), mpc_type);
                        }
                        return;
                    }

                    const ID master = model.resolve_node_id(node_2);
                    for (std::size_t i = 0; i < first->size(); ++i) {
                        mpcs.emplace_back(first->at(i), master, mpc_type);
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
