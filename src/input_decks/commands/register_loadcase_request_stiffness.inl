// register_loadcase_request_stiffness.inl — registers REQUESTSTIFFNESS within *LOADCASE

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../loadcase/linear_buckling.h"
#include "../../loadcase/linear_static.h"

namespace fem::input_decks::commands {

inline void register_loadcase_request_stiffness(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("REQUESTSTIFFNESS", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Request stiffness matrix output for supported loadcases.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("FILE").optional()
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("REQUESTSTIFFNESS must appear inside *LOADCASE");
            }

            std::string file = keys.has("FILE") ? keys.raw("FILE") : "stiffness_" + std::to_string(base->get_id()) + ".txt";

            if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                lc->stiffness_file = std::move(file);
                return;
            }
            if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                lc->stiffness_file = std::move(file);
                return;
            }

            throw std::runtime_error("REQUESTSTIFFNESS not supported for loadcase type " + parser.active_loadcase_type());
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

