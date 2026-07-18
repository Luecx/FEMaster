// register_loadcase_rebalance.inl — registers REBALANCELOADS within *LOADCASE (LinearStatic and derived)

#include <stdexcept>

#include "../parser.h"

#include "../../../loadcase/linear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_rebalance(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("REBALANCELOADS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Enable rigid-body rebalancing of external loads (sum F=M=0) for the active linear static load case.");

        // Toggle only; no additional keywords
        command.on_enter([&parser](const fem::io::dsl::Keys&) {
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

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands

