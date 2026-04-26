// register_solid_section.inl — DSL registration for *SOLIDSECTION

#include <memory>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_solid_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("SOLIDSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a material to a solid element set.");

        auto material = std::make_shared<std::string>();
        auto elset = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MATERIAL")
                    .alternative("MAT")
                    .required()
                    .doc("Material name")
                .key("ELSET")
                    .required()
                    .doc("Target element set")
                .key("ORIENTATION")
                    .optional()
                    .doc("Optional coordinate system for solid material axes")
        );

        command.on_enter([&model, material, elset, orientation](const fem::dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            model.solid_section(*elset, *material, *orientation);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
