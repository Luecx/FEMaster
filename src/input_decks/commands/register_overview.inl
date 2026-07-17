// register_overview.inl — registers OVERVIEW (prints model summary)

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model_overview.h"

namespace fem::input_decks::commands {

inline void register_overview(fem::io::dsl::Registry& registry, fem::model::Model& model) {
    registry.command("OVERVIEW", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Print a summary/overview of the current model (nodes, elements, sets, materials, sections, couplings).");

        command.on_enter([&model](const fem::io::dsl::Keys&) {
            fem::model::print_model_overview(model);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

