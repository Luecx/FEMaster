// register_loadcase_request_stgeom.inl — registers REQUESTSTGEOM for buckling loadcases

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../loadcase/linear_buckling.h"

namespace fem::input_decks::commands {

inline void register_loadcase_request_stgeom(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("REQUESTSTGEOM", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Request geometric stiffness output for LINEARBUCKLING loadcases.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("FILE").optional()
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* lc = parser.active_loadcase_as<loadcase::LinearBuckling>();
            if (!lc) {
                throw std::runtime_error("REQUESTSTGEOM only valid for LINEARBUCKLING loadcases");
            }
            lc->geom_file = keys.has("FILE") ? keys.raw("FILE") : "geom_" + std::to_string(lc->get_id()) + ".txt";
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

