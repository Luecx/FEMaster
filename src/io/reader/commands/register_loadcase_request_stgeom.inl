// register_loadcase_request_stgeom.inl — registers REQUESTSTGEOM for buckling loadcases

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../../loadcase/linear_buckling.h"

namespace fem::io::reader::commands {

inline void register_loadcase_request_stgeom(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("REQUESTSTGEOM", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Request geometric stiffness output for LINEARBUCKLING loadcases.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("FILE").optional()
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* lc = parser.active_loadcase_as<loadcase::LinearBuckling>();
            if (!lc) {
                throw std::runtime_error("REQUESTSTGEOM only valid for LINEARBUCKLING loadcases");
            }
            lc->geom_file = keys.has("FILE") ? keys.raw("FILE") : "geom_" + std::to_string(lc->get_id()) + ".txt";
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands

