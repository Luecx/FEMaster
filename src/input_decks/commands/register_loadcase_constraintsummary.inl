// register_loadcase_constraintsummary.inl — registers CONSTRAINTSUMMARY inside *LOADCASE

#include <stdexcept>

#include "../parser.h"

namespace fem::input_decks::commands {

inline void register_loadcase_constraintsummary(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("CONSTRAINTSUMMARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Enable constraint summary output for the active loadcase.");

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto* lc = parser.active_loadcase();
            if (!lc) {
                throw std::runtime_error("CONSTRAINTSUMMARY must appear inside *LOADCASE");
            }
            lc->report_constraints = true;
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

