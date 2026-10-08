/**
 * @file register_loadcase_loads.cpp
 * @brief Registers load-collector selection for active load cases.
 *
 * The LOADS child command resolves named load collectors immediately against
 * ModelData::load_cols and activates their shared condition definitions in the
 * model-owned ConditionManager. It remains available in native LOADCASE and
 * Abaqus STEP scopes independently of persistent direct condition history.
 *
 * Entering a new analysis clears only collector-derived activations. Named
 * definitions and persistent direct history remain model-owned and are unaffected.
 * Missing names are reported by the DSL before the solver runs; no analysis class
 * owns collector identifiers or performs deferred name resolution.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <string>

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

/**
 * Registers immediate named loads selection in the current model state.
 *
 * Each non-empty token must identify an existing collector. The callback expands
 * its entries directly into the current load activations in ConditionManager.
 * Definition storage and persistent direct condition history are not modified.
 *
 * @param registry Registry receiving the analysis child command.
 * @param parser Parser supplying the model and active analysis scope.
 */
void register_loadcase_loads(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("LOADS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Assign named load collectors to the active analysis step.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<std::string, 16>().name("LOAD").desc("Load collector names")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&parser](const std::array<std::string, 16>& names) {
                    auto* base = parser.active_loadcase();
                    logging::error(base != nullptr,
                        "LOADS must appear inside an active analysis step");

                    // Resolve every user-supplied name now. The solver later reads
                    // only the resulting ConditionManager state.
                    auto& model_data = *parser.model()._data;
                    for (const auto& name : names) {
                        if (name.empty()) {
                            continue;
                        }

                        logging::error(model_data.load_cols.has(name),
                            "LOADS: collector ", name, " does not exist");
                        const auto collector = model_data.load_cols.get(name);
                        logging::error(collector != nullptr,
                            "LOADS: collector ", name, " is not initialized");

                        // Select reusable definitions for this analysis only.
                        // The parser suppresses repeated (family, pointer)
                        // selections before the temporary manager insertion.
                        for (const auto& condition : *collector) {
                            parser.select_collector_condition(bc::CLOAD, condition);
                        }
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
