// register_material.inl — DSL registration for *MATERIAL

#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::input_decks::commands {

inline void register_material(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("MATERIAL", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Activate or create a material definition.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MATERIAL")
                    .alternative("NAME")
                    .required()
                    .doc("Identifier of the material entry")
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string& name = keys.raw("MATERIAL");
            model._data->materials.activate(name);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
