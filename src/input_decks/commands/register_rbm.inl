// register_rbm.inl — registers *RBM

#include <string>
#include <stdexcept>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_rbm(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("RBM", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Add rigid-body-motion suppression equations for an element set. "
            "ELSET defaults to EALL."
        );

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("ELSET")
                    .alternative("SET")
                    .optional("EALL")
                    .doc("Target element set used to build RBM equations")
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            if (keys.has("MAX_POINTS")) {
                throw std::runtime_error("RBM key 'MAX_POINTS' is no longer supported");
            }
            const std::string elset = keys.raw("ELSET");
            model.add_rbm(elset);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
