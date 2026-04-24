// register_shell_section.inl — DSL registration for *SHELLSECTION

#include <memory>
#include <string>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_shell_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("SHELLSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a shell section thickness to an element set.");

        // Persistent per-command state
        auto material = std::make_shared<std::string>();
        auto elset    = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required().doc("Material name")
                .key("ELSET").required().doc("Target element set")
                .key("ORIENTATION").optional().doc("Optional coordinate system for shell n1/n2 material/resultant axes")
        );

        // Capture shared_ptrs BY VALUE so they're valid later
        command.on_enter([material, elset, orientation](const fem::dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset    = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::Precision>().name("THICKNESS").desc("Shell thickness")
                        .on_missing(fem::Precision{1}).on_empty(fem::Precision{1})
                )
                .bind([&model, material, elset, orientation](fem::Precision thickness) {
                    model.shell_section(*elset, *material, thickness, *orientation);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
