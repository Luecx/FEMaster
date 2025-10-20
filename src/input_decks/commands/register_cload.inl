// register_cload.inl — DSL registration for *CLOAD
#pragma once
/**
 * @file register_cload.inl
 * @brief Register *CLOAD command.
 *
 * Concentrated nodal loads: forces (Fx,Fy,Fz) and moments (Mx,My,Mz)
 * applied in the global coordinate system about the node.
 */

#include <array>
#include <stdexcept>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_cload(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("CLOAD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));

        command.doc(
            "Applies concentrated nodal loads in the global CS at the node location. "
            "Forces (Fx,Fy,Fz) and moments (Mx,My,Mz) follow the right-handed sign convention. "
            "Multiple lines targeting the same node/nodes are accumulated."
        );


        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").required()
                    .doc("Target load collector to which these nodal loads are added.")
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& collector = keys.raw("LOAD_COLLECTOR");
            model._data->load_cols.activate(collector);
        });

        command.variant(
            fem::dsl::Variant::make()
                .doc("Each data line defines a target and six components: Fx, Fy, Fz, Mx, My, Mz.")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1)) // one or more lines
                        .pattern(
                            fem::dsl::Pattern::make()
                                .one<std::string>()
                                    .name("TARGET")
                                    .desc("Node set name or single node id (integer).")
                                .fixed<fem::Precision, 6>()
                                    .name("LOAD")
                                    .desc("Fx, Fy, Fz, Mx, My, Mz (global).")
                                    .on_missing(fem::Precision{0})
                                    .on_empty  (fem::Precision{0})
                        )
                        .bind([&model](const std::string& target,
                                       const std::array<fem::Precision, 6>& values) {
                            fem::Vec6 load;
                            load << values[0], values[1], values[2],
                                    values[3], values[4], values[5];

                            // node set target?
                            if (model._data->node_sets.has(target)) {
                                model.add_cload(target, load);
                                return;
                            }

                            // single node id?
                            try {
                                const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                                model.add_cload(id, load);
                            } catch (const std::exception&) {
                                throw std::runtime_error(
                                    "CLOAD target '" + target +
                                    "' is neither an existing node set nor a valid node id."
                                );
                            }
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
