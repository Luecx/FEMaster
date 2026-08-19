/**
 * @file register_tie.inl
 * @brief Registers tie constraints.
 *
 * The parser resolves whether the slave is nodal or surface-based and whether
 * the master uses surfaces or lines. It then constructs the exact Tie variant
 * and stores it directly in ModelData. No Model-side constraint factory or type
 * dispatcher participates in this path.
 *
 * Search distance and optional initial adjustment are retained by the selected
 * tie formulation. Projection, interpolation and constraint-equation assembly
 * use the compiled geometry later when a load case collects constraints.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../../constraints/types/tie.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_tie(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("TIE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER").required()
                .key("SLAVE").required()
                .key("ADJUST").optional("NO").allowed({"NO", "YES"})
                .key("DISTANCE").required()
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string master = keys.raw("MASTER");
            const std::string slave  = keys.raw("SLAVE");
            const Precision distance = keys.get<Precision>("DISTANCE");
            const bool adjust        = keys.get<bool>("ADJUST");

            logging::error(model._data->compiled,
                "TIE: constraints require a compiled model");

            const auto slave_nodes = model._data->node_sets.has(slave)
                ? model._data->node_sets.get(slave)
                : nullptr;
            const auto slave_surfaces = model._data->surface_sets.has(slave)
                ? model._data->surface_sets.get(slave)
                : nullptr;

            logging::error((slave_nodes && slave_nodes->size() > 0)
                        || (slave_surfaces && slave_surfaces->size() > 0),
                "TIE: slave set ", slave, " is not a non-empty node or surface set");

            const auto master_surfaces = model._data->surface_sets.has(master)
                ? model._data->surface_sets.get(master)
                : nullptr;
            const auto master_lines = model._data->line_sets.has(master)
                ? model._data->line_sets.get(master)
                : nullptr;

            logging::error((master_surfaces && master_surfaces->size() > 0)
                        || (master_lines && master_lines->size() > 0),
                "TIE: master set ", master, " is not a non-empty surface or line set");

            if (master_surfaces && master_surfaces->size() > 0) {
                if (slave_nodes && slave_nodes->size() > 0) {
                    model._data->ties.emplace_back(master_surfaces, slave_nodes, distance, adjust);
                } else {
                    model._data->ties.emplace_back(master_surfaces, slave_surfaces, distance, adjust);
                }
            } else if (slave_nodes && slave_nodes->size() > 0) {
                model._data->ties.emplace_back(master_lines, slave_nodes, distance, adjust);
            } else {
                model._data->ties.emplace_back(master_lines, slave_surfaces, distance, adjust);
            }
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
