// register_heading.inl — DSL registration for *HEADING

#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::input_decks::commands {

inline void register_heading(fem::io::dsl::Registry& registry) {
    registry.command("HEADING", [](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Accept and ignore optional model heading text.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TEXT").desc("Ignored heading line")
                )
                .bind([](const std::string&) {})
            )
        );
    });
}

} // namespace fem::input_decks::commands
