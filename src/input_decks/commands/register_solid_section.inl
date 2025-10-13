// register_solid_section.inl — DSL registration for *SOLIDSECTION

#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_solid_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("SOLIDSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a material to a solid element set.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MATERIAL")
                    .alternative("MAT")
                    .required()
                    .doc("Material name")
                .key("ELSET")
                    .required()
                    .doc("Target element set")
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& material = keys.raw("MATERIAL");
            const std::string& elset    = keys.raw("ELSET");
            model.solid_section(elset, material);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
