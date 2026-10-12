/**
 * @file register_loadcase_collectors.cpp
 * @brief Registers named load and support collector selection for analyses.
 *
 * LOADS and SUPPORTS share the same DSL structure and selection semantics.
 * Named definitions are resolved against their respective model registries
 * and selected for the current analysis without modifying direct condition
 * history. The parser inserts selected definitions only during the solve.
 */

#include "register_functions.h"
#include "../../dsl/registry.h"
#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../model/model.h"

#include <array>
#include <string>

namespace fem::io::reader::commands {
namespace {

void register_collector_selection(dsl::Registry& registry, Parser& parser,
                                  bc::ConditionFamily family) {
    const bool supports = family == bc::SUPPORT;
    const std::string command_name = supports ? "SUPPORTS" : "LOADS";

    registry.command(command_name, [&parser, family, supports, command_name](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({
            "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"
        }));
        command.doc(supports
            ? "Assign named support collectors to the active analysis step."
            : "Assign named load collectors to the active analysis step.");

        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .fixed<std::string, 16>()
                        .name(supports ? "SUPP" : "LOAD")
                        .desc(supports ? "Support collector names" : "Load collector names")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&parser, family, supports, command_name](const std::array<std::string, 16>& names) {
                    logging::error(parser.active_loadcase() != nullptr,
                        command_name, " must appear inside an active analysis step");

                    auto& model_data = *parser.model()._data;
                    auto& collectors = supports ? model_data.supp_cols : model_data.load_cols;

                    for (const auto& name : names) {
                        if (name.empty()) continue;

                        logging::error(collectors.has(name),
                            command_name, ": collector ", name, " does not exist");
                        const auto collector = collectors.get(name);
                        logging::error(collector != nullptr,
                            command_name, ": collector ", name, " is not initialized");

                        for (const auto& condition : *collector) {
                            parser.select_collector_condition(family, condition);
                        }
                    }
                })
            )
        );
    });
}

} // namespace

void register_loadcase_loads(dsl::Registry& registry, Parser& parser) {
    register_collector_selection(registry, parser, bc::CLOAD);
}

void register_loadcase_supports(dsl::Registry& registry, Parser& parser) {
    register_collector_selection(registry, parser, bc::SUPPORT);
}

} // namespace fem::io::reader::commands
