// register_loadcase_constraintsummary.inl — registers CONSTRAINTSUMMARY inside *LOADCASE

#include "../parser.h"

#include "../../../core/logging.h"

namespace fem::io::reader::commands {

inline void register_loadcase_constraintsummary(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("CONSTRAINTSUMMARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Enable constraint summary output for the active loadcase.");

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto* lc = parser.active_loadcase();
            logging::error(lc != nullptr,
                "CONSTRAINTSUMMARY must appear inside *LOADCASE");
            lc->report_constraints = true;
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands

