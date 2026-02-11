// register_loadcase_rebalance.inl — registers REBALANCELOADS within *LOADCASE (LinearStatic and derived)

#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/linear_static.h"

namespace fem::input_decks::commands {

inline void register_loadcase_rebalance(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("REBALANCELOADS", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Enable rigid-body rebalancing of external loads (sum F=M=0) for the active linear static load case.");

        // Toggle only; no additional keywords
        command.on_enter([&parser](const fem::dsl::Keys&) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("REBALANCELOADS must appear inside *LOADCASE");
            }
            if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                lc->rebalance_loads = true;
                return;
            }
            throw std::runtime_error("REBALANCELOADS is only supported for linear static load cases");
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

