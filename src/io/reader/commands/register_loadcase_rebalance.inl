/**
 * @file register_loadcase_rebalance.inl
 * @brief Registers rigid-body load rebalancing for linear static analyses.
 *
 * The flag-style `REBALANCELOADS` command enables compensation of the active
 * external load system so its resultant force and moment vanish. It applies to
 * `LinearStatic` and derived formulations and only changes the corresponding
 * analysis setting.
 *
 * Evaluation of resultants, selection of compensating nodal loads and stiffness
 * solution remain part of the linear static implementation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_rebalance(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("REBALANCELOADS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Enable rigid-body rebalancing of external loads (sum F=M=0) for the active linear static load case.");

        // Toggle only; no additional keywords
        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto* base = parser.active_loadcase();
            logging::error(base != nullptr,
                "REBALANCELOADS must appear inside *LOADCASE");

            auto* lc = dynamic_cast<loadcase::LinearStatic*>(base);
            logging::error(lc != nullptr,
                "REBALANCELOADS is only supported for linear static load cases");
            lc->rebalance_loads = true;
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
