/**
 * @file register_tie.cpp
 * @brief Registers tie constraints.
 *
 * The parser resolves whether the slave is nodal or surface-based and whether
 * the master uses surfaces or lines. It then constructs the exact Tie variant
 * and stores it directly in ModelData. No Model-side constraint factory or type
 * dispatcher participates in this path.
 *
 * `TIE, ADJUST` is applied immediately after construction against the compiled
 * reference geometry. Constraint-equation assembly remains deferred until a
 * load case collects constraints.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"

#include <memory>

#include "../../dsl/registry.h"

#include "../../../constraints/types/tie.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

void register_tie(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("TIE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional()
                .key("MASTER").optional()
                .key("SLAVE").optional()
                .key("ADJUST").optional("NO").allowed({"NO", "YES"})
                .key("DISTANCE").alternative("POSITIONTOLERANCE").required()
        );
        const auto distance = std::make_shared<Precision>(0);
        const auto adjust   = std::make_shared<bool>(false);

        const auto add_tie = [&model](const std::string& slave, const std::string& master,
                                      Precision max_distance, bool do_adjust) {
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
                    model._data->ties.emplace_back(master_surfaces, slave_nodes, max_distance, do_adjust);
                } else {
                    model._data->ties.emplace_back(master_surfaces, slave_surfaces, max_distance, do_adjust);
                }
            } else if (slave_nodes && slave_nodes->size() > 0) {
                model._data->ties.emplace_back(master_lines, slave_nodes, max_distance, do_adjust);
            } else {
                model._data->ties.emplace_back(master_lines, slave_surfaces, max_distance, do_adjust);
            }

            model._data->ties.back().adjust_geometry(*model._data);
        };

        command.on_enter([&model, distance, adjust, add_tie](const fem::io::dsl::Keys& keys) {
            logging::error(model._data->compiled,
                "TIE: constraints require a compiled model");
            logging::error(keys.has("MASTER") == keys.has("SLAVE"),
                "TIE: MASTER and SLAVE must either both be specified or both omitted");

            *distance = keys.get<Precision>("DISTANCE");
            *adjust   = keys.get<bool>("ADJUST");
            if (keys.has("MASTER")) {
                add_tie(keys.raw("SLAVE"), keys.raw("MASTER"), *distance, *adjust);
            }
        });

        // Abaqus format: one or more data lines, each containing slave, master.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::negate(
                fem::io::dsl::Condition::key_present("MASTER")))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("SLAVE")
                    .one<std::string>().name("MASTER")
                )
                .bind([add_tie, distance, adjust](const std::string& slave,
                                                  const std::string& master) {
                    add_tie(slave, master, *distance, *adjust);
                })
            )
        );

        // Preserve existing FEMaster decks that supply both names on the keyword.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_present("MASTER"))
        );
    });
}

} // namespace fem::io::reader::commands
