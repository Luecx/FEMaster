// register_loadcase_request_stiffness.inl — registers REQUESTSTIFFNESS within *LOADCASE

#include <string>

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_request_stiffness(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("REQUESTSTIFFNESS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Request stiffness matrix output for supported loadcases.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("FILE").optional()
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            logging::error(base != nullptr,
                "REQUESTSTIFFNESS must appear inside *LOADCASE");

            std::string file = keys.has("FILE") ? keys.raw("FILE") : "stiffness_" + std::to_string(base->get_id()) + ".txt";

            if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                lc->stiffness_file = std::move(file);
                return;
            }
            if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                lc->stiffness_file = std::move(file);
                return;
            }
            if (auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base)) {
                lc->stiffness_file = std::move(file);
                return;
            }

            logging::error(false,
                "REQUESTSTIFFNESS not supported for loadcase type ", parser.active_loadcase_type());
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
