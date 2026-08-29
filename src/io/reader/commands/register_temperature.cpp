/**
 * @file register_temperature.cpp
 * @brief Registers prescribed nodal temperatures as support conditions.
 *
 * The root- and assembly-level `TEMPERATURE` command targets compiled nodes or
 * node regions and stores each prescription in the named common support
 * collector. Mechanical and thermal support definitions share collector
 * ownership through `bc::SupportInterface`; thermal load cases later select
 * only `bc::Temperature` entries.
 *
 * Temperature values are absolute scalar primary variables. They are therefore
 * essential boundary conditions rather than thermal loads and do not contribute
 * to a heat-source right-hand side.
 *
 * @see bc::Temperature
 * @see bc::SupportCollector
 * @see model::Model::add_support
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <memory>
#include <string>
#include <utility>

#include "../../../bc/thermal_temp.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers prescribed absolute temperatures for compiled node targets.
 *
 * Command entry activates the requested shared support collector. Every data row
 * resolves its target as an existing compiled node region or as one scalar node
 * reference and transfers a polymorphic `bc::Temperature` to that collector.
 *
 * @param registry DSL registry receiving the temperature command.
 * @param model Compiled model providing node resolution and support storage.
 */
void register_temperature(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("TEMPERATURE", [&](fem::io::dsl::Command& command) {
        // Prescribed temperatures operate only on compiled assembly nodes
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Prescribe absolute nodal temperatures in a named support collector.");

        // Select the common collector receiving all rows of this command
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SUPPORT_COLLECTOR")
                    .alternative("SUPPORTCOLLECTOR")
                    .alternative("NAME")
                    .required()
                    .doc("Common support collector receiving the temperatures")
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->supp_cols.activate(keys.raw("SUPPORT_COLLECTOR"));
        });

        // Resolve node-set or scalar-node targets and store absolute values
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled node or node-set reference")
                    .one<fem::Precision>().name("TEMPERATURE").desc("Absolute prescribed temperature")
                )
                .bind([&model](const std::string& target, fem::Precision temperature) {
                    model::NodeRegion::Ptr region;
                    if (model._data->node_sets.has(target)) {
                        region = model._data->node_sets.get(target);
                    } else {
                        region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(model.compiled_node_id(target));
                    }

                    model.add_support(std::make_shared<bc::Temperature>(
                        std::move(region),
                        temperature
                    ));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
