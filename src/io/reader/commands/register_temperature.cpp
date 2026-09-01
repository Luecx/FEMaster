/**
 * @file register_temperature.cpp
 * @brief Registers prescribed-temperature Dirichlet boundary conditions.
 *
 * `TEMPERATURE` is the thermal essential-boundary-condition counterpart of
 * structural `SUPPORT`. Each data row assigns one absolute scalar temperature
 * to a compiled node or node set and stores the resulting `bc::Temperature` in
 * the named common support collector. Structural and thermal Dirichlet
 * conditions therefore share collector ownership while individual analyses can
 * select only the concrete condition types relevant to their field.
 *
 * Native syntax:
 *
 * @code
 * *TEMPERATURE, SUPPORT_COLLECTOR=THERMAL_BCS
 * NSET_HOT, 373.15
 * 42,       293.15
 * @endcode
 *
 * Temperature is a primary unknown and consequently contributes constraint
 * equations rather than entries to the thermal right-hand side.
 *
 * @see bc::Temperature
 * @see bc::Dirichlet
 * @see bc::SupportCollector
 * @see model::Model::add_support
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <memory>
#include <string>
#include <utility>

#include "../../../bc/dirichlet/temperature.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * @brief Registers absolute nodal temperatures in the native DSL.
 *
 * The command is valid after assembly compilation at root or assembly scope.
 * Entering the command activates the requested support collector. Each row then
 * resolves `TARGET` first as a compiled node set and otherwise as one compiled
 * node reference. The parsed scalar is passed unchanged to `bc::Temperature`.
 *
 * @param registry DSL registry receiving the `TEMPERATURE` command.
 * @param model Compiled model used for target resolution and support storage.
 */
void register_temperature(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("TEMPERATURE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Prescribe absolute nodal temperature in a named support collector.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SUPPORT_COLLECTOR")
                    .alternative("SUPPORT COLLECTOR")
                    .alternative("NAME")
                    .required()
                    .doc("Support collector receiving the prescribed temperatures")
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->supp_cols.activate(keys.raw("SUPPORT_COLLECTOR"));
        });

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
                        std::move(region), temperature));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
