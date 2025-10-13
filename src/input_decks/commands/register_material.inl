// register_material.inl — DSL registration for *MATERIAL

#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::input_decks::commands {

inline void register_material(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("MATERIAL", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Activate or create a material definition.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MATERIAL")
                    .alternative("NAME")
                    .required()
                    .doc("Identifier of the material entry")
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& name = keys.raw("MATERIAL");
            model._data->materials.activate(name);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
