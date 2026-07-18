// register_heading.inl — DSL registration for *HEADING

#include <array>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_heading(fem::io::dsl::Registry& registry) {
    registry.command("HEADING", [](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Accept and ignore optional model heading text.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<std::string, 64>().name("TEXT").desc("Ignored heading line").on_missing(std::string{"_"})
                )
                .bind([](const std::array<std::string, 64>&) {})
            )
        );
    });
}

} // namespace fem::io::reader::commands
